#!/usr/bin/env python3
"""Split a registered DBiT image and segment cells in every spatial spot."""

from __future__ import annotations

import argparse
import io
import sys
from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile
from skimage.segmentation import mark_boundaries
from stardist.models import StarDist2D
from tqdm import tqdm


ORIENTATION_CHOICES = ("normal", "horizontal", "vertical", "rotate")
IMAGE_SUFFIXES = (".tif", ".tiff", ".png", ".jpg")

BACKGROUND_THRESHOLD = 8
LOCAL_DENSITY_KERNEL = 201
MIN_LOCAL_SIGNAL_FRACTION = 0.02
MORPHOLOGY_KERNEL = 11
MIN_COMPONENT_FRACTION = 0.0001
MIN_COMPONENT_RELATIVE_TO_LARGEST = 0.001
MIN_TISSUE_FRACTION = 0.03


@dataclass(frozen=True)
class GridConfig:
    x_spots_number: int
    y_spots_number: int
    length_spot: int
    interval: int
    pixel_length: float


@dataclass(frozen=True)
class StarDistConfig:
    top_value: int
    number_of_top_values: int
    prob_thresh: float
    nms_thresh: float
    model_name: str


def str_to_bool(value):
    """Convert a common command-line Boolean string to bool."""
    if isinstance(value, bool):
        return value
    normalized = value.lower()
    if normalized in ("yes", "true", "t", "y", "1"):
        return True
    if normalized in ("no", "false", "f", "n", "0"):
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got: {value}")


def normalize_orientation(value: str) -> str:
    orientation = value.lower()
    if orientation not in ORIENTATION_CHOICES:
        raise argparse.ArgumentTypeError(
            "orientation must be one of "
            f"{', '.join(ORIENTATION_CHOICES)}, got: {value}"
        )
    return orientation


def _grayscale(image: np.ndarray) -> np.ndarray:
    if image.ndim == 2:
        return image
    if image.ndim != 3:
        raise ValueError(f"Unsupported image shape: {image.shape}")
    if image.shape[2] == 1:
        return image[:, :, 0]
    return cv2.cvtColor(image[:, :, :3], cv2.COLOR_BGR2GRAY)


def _to_uint8(image: np.ndarray) -> np.ndarray:
    if image.dtype == np.uint8:
        return image
    normalized = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX)
    return normalized.astype(np.uint8)


def _remove_small_regions(mask: np.ndarray) -> np.ndarray:
    contours, _ = cv2.findContours(
        mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
    )
    if not contours:
        return np.zeros_like(mask)
    areas = [cv2.contourArea(contour) for contour in contours]
    min_area = max(
        64,
        mask.size * MIN_COMPONENT_FRACTION,
        max(areas) * MIN_COMPONENT_RELATIVE_TO_LARGEST,
    )
    cleaned = np.zeros_like(mask)
    retained = [
        contour for contour, area in zip(contours, areas) if area >= min_area
    ]
    cv2.drawContours(cleaned, retained, -1, 255, thickness=cv2.FILLED)
    return cleaned


def generate_tissue_mask(image_path: Path, output_path: Path) -> np.ndarray:
    """Generate and save a coarse whole-image tissue mask."""
    image = cv2.imread(str(image_path), cv2.IMREAD_UNCHANGED)
    if image is None:
        raise ValueError(f"Unable to read image: {image_path}")

    intensity = _to_uint8(_grayscale(image))
    signal = np.where(intensity > BACKGROUND_THRESHOLD, 255, 0).astype(np.uint8)
    local_density = cv2.boxFilter(
        signal,
        ddepth=-1,
        ksize=(LOCAL_DENSITY_KERNEL, LOCAL_DENSITY_KERNEL),
        normalize=True,
    )
    density_threshold = round(255 * MIN_LOCAL_SIGNAL_FRACTION)
    mask = np.where(local_density >= density_threshold, 255, 0).astype(np.uint8)
    kernel = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (MORPHOLOGY_KERNEL, MORPHOLOGY_KERNEL)
    )
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
    mask = _remove_small_regions(mask)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(output_path), mask):
        raise RuntimeError(f"Failed to write tissue mask: {output_path}")
    return mask


def split_image(
    image_path: Path,
    grid: GridConfig,
    split_path: Path,
    result_path: Path,
    put_text: bool,
    font_size: int,
    orientation: str,
    swap_xy: bool,
    tissue_mask: np.ndarray,
    write_tiles: bool,
) -> dict[tuple[int, int], bool]:
    """Classify every grid spot and optionally write its cropped image tile."""
    image = cv2.imread(str(image_path))
    if image is None:
        raise ValueError(f"Unable to read image: {image_path}")

    x_starts = [
        int(x * (grid.length_spot + grid.interval) / grid.pixel_length)
        for x in range(grid.x_spots_number)
    ]
    y_starts = [
        int(y * (grid.length_spot + grid.interval) / grid.pixel_length)
        for y in range(grid.y_spots_number)
    ]
    spot_pixels = int(grid.length_spot / grid.pixel_length)
    x_ends = [start + spot_pixels for start in x_starts]
    y_ends = [start + spot_pixels for start in y_starts]

    required_width = max(x_ends)
    required_height = max(y_ends)
    image_height, image_width = image.shape[:2]
    if image_width < required_width or image_height < required_height:
        raise ValueError(
            "Image is smaller than the configured DBiT grid: "
            f"image={image_width}x{image_height}, "
            f"required={required_width}x{required_height}"
        )
    if tissue_mask.shape[:2] != image.shape[:2]:
        raise ValueError(
            f"Tissue mask shape {tissue_mask.shape[:2]} does not match "
            f"image shape {image.shape[:2]}"
        )

    output_x_count = grid.y_spots_number if swap_xy else grid.x_spots_number
    output_y_count = grid.x_spots_number if swap_xy else grid.y_spots_number
    reverse_x = orientation in ("horizontal", "rotate")
    reverse_y = orientation in ("vertical", "rotate")
    tissue_spots: dict[tuple[int, int], bool] = {}

    description = (
        "Splitting image and classifying tissue"
        if write_tiles
        else "Classifying tissue"
    )
    for x in tqdm(range(output_x_count), desc=description):
        for y in range(output_y_count):
            mapped_x = output_x_count - 1 - x if reverse_x else x
            mapped_y = output_y_count - 1 - y if reverse_y else y
            crop_x, crop_y = (
                (mapped_y, mapped_x) if swap_xy else (mapped_x, mapped_y)
            )
            x_start, x_end = x_starts[crop_x], x_ends[crop_x]
            y_start, y_end = y_starts[crop_y], y_ends[crop_y]

            mask_crop = tissue_mask[y_start:y_end, x_start:x_end]
            tissue_fraction = float((mask_crop > 0).sum()) / float(mask_crop.size)
            tissue_spots[(x, y)] = tissue_fraction >= MIN_TISSUE_FRACTION

            if write_tiles:
                tile = image[y_start:y_end, x_start:x_end]
                tile_path = split_path / f"{x}_{y}.tif"
                if not cv2.imwrite(str(tile_path), tile):
                    raise RuntimeError(f"Failed to write image tile: {tile_path}")
                cv2.rectangle(
                    image, (x_start, y_start), (x_end, y_end), (255, 0, 0), 2
                )
                if put_text:
                    cv2.putText(
                        image,
                        f"{x}_{y}",
                        (x_start, y_end),
                        cv2.FONT_HERSHEY_SIMPLEX,
                        font_size,
                        (0, 0, 255),
                        1,
                    )

    if write_tiles:
        overview_path = result_path / "result.png"
        if not cv2.imwrite(str(overview_path), image):
            raise RuntimeError(f"Failed to write grid overview: {overview_path}")
    return tissue_spots


_MODEL_CACHE: dict[str, StarDist2D] = {}


def _get_model(model_name: str) -> StarDist2D:
    if model_name not in _MODEL_CACHE:
        stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            _MODEL_CACHE[model_name] = StarDist2D.from_pretrained(model_name)
        finally:
            sys.stdout = stdout
    return _MODEL_CACHE[model_name]


def predict_tile(
    image_path: Path,
    mask_path: Path,
    label_path: Path,
    config: StarDistConfig,
) -> np.ndarray | None:
    """Run StarDist for one tile and write its binary mask and boundary overlay."""
    original = tifffile.imread(image_path)
    image = original.copy()
    intensity = image[:, :, 1] if image.ndim == 3 and image.shape[-1] >= 2 else image
    flat_intensity = intensity.ravel()
    top_count = min(config.number_of_top_values, flat_intensity.size)
    top_values = np.partition(flat_intensity, -top_count)[-top_count:]
    if np.mean(top_values) <= config.top_value:
        return None

    if image.ndim == 3 and image.shape[-1] >= 2:
        normalized = image[:, :, 1]
    elif image.ndim == 3 and image.shape[-1] == 1:
        normalized = image[:, :, 0]
    else:
        normalized = image
    normalized = normalized.astype(np.float32)
    if normalized.max() > 0:
        normalized /= normalized.max()

    model = _get_model(config.model_name)
    labels, _ = model.predict_instances(
        normalized,
        prob_thresh=config.prob_thresh,
        nms_thresh=config.nms_thresh,
    )

    binary_mask = (labels > 0).astype(np.uint8) * 255
    tifffile.imwrite(mask_path, binary_mask)

    if original.ndim == 3 and original.shape[-1] == 3:
        visualization = original.astype(np.float32)
    elif original.ndim == 3 and original.shape[-1] == 1:
        green = original[:, :, 0].astype(np.float32)
        zeros = np.zeros_like(green)
        visualization = np.stack([zeros, green, zeros], axis=2)
    else:
        grayscale = original.astype(np.float32)
        visualization = np.repeat(grayscale[:, :, np.newaxis], 3, axis=2)
    if visualization.max() > 0:
        visualization /= visualization.max()

    overlay = mark_boundaries(
        visualization, labels, color=(1, 0, 0), mode="thick"
    )
    tifffile.imwrite(label_path, (overlay * 255).astype(np.uint8))
    return labels


def run_segmentation(
    image_path: Path,
    result_path: Path,
    grid: GridConfig,
    stardist: StarDistConfig,
    put_text: bool,
    font_size: int,
    orientation: str,
    swap_xy: bool,
) -> pd.DataFrame:
    """Split the registered image, predict cells, and write cell_num_area.csv."""
    split_path = result_path / "split"
    mask_path = result_path / "mask"
    label_path = result_path / "label"
    for directory in (split_path, mask_path, label_path):
        directory.mkdir(parents=True, exist_ok=True)

    tissue_mask = generate_tissue_mask(image_path, result_path / "tissue_mask.png")
    split_files = [
        path for path in split_path.iterdir() if path.suffix.lower() in IMAGE_SUFFIXES
    ]
    orientation_file = split_path / ".orientation"
    transform_key = f"orientation={orientation};swap_xy={swap_xy}"
    previous_orientation = (
        orientation_file.read_text(encoding="utf-8").strip()
        if orientation_file.exists()
        else ""
    )
    expected_split_count = grid.x_spots_number * grid.y_spots_number
    write_tiles = not (
        len(split_files) == expected_split_count
        and previous_orientation == transform_key
    )
    if not write_tiles:
        print(
            f"✓ Split images already exist in {split_path} with "
            f"{transform_key}, skipping splitting step."
        )

    tissue_spots = split_image(
        image_path=image_path,
        grid=grid,
        split_path=split_path,
        result_path=result_path,
        put_text=put_text,
        font_size=font_size,
        orientation=orientation,
        swap_xy=swap_xy,
        tissue_mask=tissue_mask,
        write_tiles=write_tiles,
    )
    print(f"Tissue spots: {sum(tissue_spots.values())}/{len(tissue_spots)}")
    if write_tiles:
        orientation_file.write_text(transform_key, encoding="utf-8")

    records = []
    for tile_path in sorted(split_path.iterdir()):
        if tile_path.suffix.lower() not in IMAGE_SUFFIXES:
            continue
        x, y = (int(value) for value in tile_path.stem.split("_"))
        try:
            labels = predict_tile(
                tile_path,
                mask_path / tile_path.name,
                label_path / tile_path.name,
                stardist,
            )
            if labels is None:
                record = {
                    "x": x,
                    "y": y,
                    "num_cells": 0,
                    "area": [],
                    "status": "skipped",
                    "in_tissue": tissue_spots[(x, y)],
                }
                print(f"✗ Skipped: {tile_path.name} - no cells")
            else:
                num_cells = int(np.max(labels))
                areas = [
                    int(np.sum(labels == cell_id))
                    for cell_id in range(1, num_cells + 1)
                ]
                record = {
                    "x": x,
                    "y": y,
                    "num_cells": num_cells,
                    "area": areas,
                    "status": "predicted",
                    "in_tissue": tissue_spots[(x, y)],
                }
                print(f"✓ Processed: {tile_path.name} - {num_cells} cells")
        except Exception as error:
            record = {
                "x": x,
                "y": y,
                "num_cells": 0,
                "area": [],
                "status": f"error: {error}",
                "in_tissue": tissue_spots[(x, y)],
            }
            print(f"✗ Error processing {tile_path.name}: {error}")
        records.append(record)

    result = pd.DataFrame.from_records(
        records,
        columns=["x", "y", "num_cells", "area", "status", "in_tissue"],
    )
    result.to_csv(result_path / "cell_num_area.csv", index=False)
    print("\n=== Processing Complete ===")
    print(f"Total processed: {len(result)}")
    print(f'Successfully predicted: {(result["status"] == "predicted").sum()}')
    print(f'Skipped/Failed: {(result["status"] != "predicted").sum()}')
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split a registered DBiT image and segment cells with StarDist."
    )
    parser.add_argument("-ip", "--image_path", required=True)
    parser.add_argument("-r", "--result_path", required=True)
    parser.add_argument("-x", "--x_spots_number", type=int, default=50)
    parser.add_argument("-y", "--y_spots_number", type=int, default=50)
    parser.add_argument("-l", "--length_spot", type=int, default=50)
    parser.add_argument("-i", "--interval", type=int, default=50)
    parser.add_argument("-p", "--pixel_length", type=float, default=0.294)
    parser.add_argument("-t", "--put_text", type=str_to_bool, default=True)
    parser.add_argument("-fs", "--font_size", type=int, default=1)
    parser.add_argument("-top_value", "--top_value", type=int, default=50)
    parser.add_argument(
        "-number_of_top_values", "--number_of_top_values", type=int, default=1500
    )
    parser.add_argument("-m", "--model_name", default="2D_versatile_fluo")
    parser.add_argument("-pt", "--prob_thresh", type=float, default=0.5)
    parser.add_argument("-nt", "--nms_thresh", type=float, default=0.6)
    parser.add_argument(
        "--orientation", type=normalize_orientation, default="normal"
    )
    parser.add_argument("--swap_xy", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    grid = GridConfig(
        args.x_spots_number,
        args.y_spots_number,
        args.length_spot,
        args.interval,
        args.pixel_length,
    )
    stardist = StarDistConfig(
        args.top_value,
        args.number_of_top_values,
        args.prob_thresh,
        args.nms_thresh,
        args.model_name,
    )
    run_segmentation(
        image_path=Path(args.image_path),
        result_path=Path(args.result_path),
        grid=grid,
        stardist=stardist,
        put_text=args.put_text,
        font_size=args.font_size,
        orientation=args.orientation,
        swap_xy=args.swap_xy,
    )


if __name__ == "__main__":
    main()
