"""Composite spatial frames onto the frame region of a full-resolution image."""

import argparse
from pathlib import Path

from PIL import Image


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge spatial frames with the masked region of an image."
    )
    parser.add_argument("--image", required=True, type=Path)
    parser.add_argument("--mask", required=True, type=Path)
    parser.add_argument("--x_spots_number", required=True, type=int)
    parser.add_argument("--y_spots_number", required=True, type=int)
    parser.add_argument("--length_spot", required=True, type=int)
    parser.add_argument("--interval", required=True, type=int)
    parser.add_argument("--pixel_length", required=True, type=float)
    parser.add_argument(
        "--frame",
        action="append",
        required=True,
        type=Path,
        help="Exact path to an overlay image; may be specified more than once",
    )
    return parser.parse_args()


Image.MAX_IMAGE_PIXELS = None


def frame_bbox(mask_path: Path) -> tuple[int, int, int, int]:
    with Image.open(mask_path) as mask:
        if "A" in mask.getbands():
            region = mask.getchannel("A")
            bbox = region.getbbox()
            if bbox is not None and region.getextrema()[0] == 0:
                return bbox
        bbox = mask.convert("L").getbbox()
    if bbox is None:
        raise ValueError(f"Frame mask contains no nonzero region: {mask_path}")
    return bbox


def clipped_axis_window(
    start: int,
    extent: int,
    canvas_extent: int,
    expected_extent: int,
    axis_name: str,
) -> tuple[int, int]:
    if extent >= expected_extent:
        return 0, expected_extent

    touches_low = start == 0
    touches_high = start + extent == canvas_extent
    if touches_low == touches_high:
        raise ValueError(
            f"Visible frame {axis_name} extent {extent} is smaller than the "
            f"configured extent {expected_extent}, but it does not touch exactly "
            f"one image boundary; the clipped-frame position is ambiguous"
        )
    if touches_low:
        return expected_extent - extent, expected_extent
    return 0, extent


def frame_crop_box(
    bbox: tuple[int, int, int, int],
    image_size: tuple[int, int],
    expected_size: tuple[int, int],
    frame_size: tuple[int, int],
) -> tuple[int, int, int, int]:
    left, top, right, bottom = bbox
    visible_width = right - left
    visible_height = bottom - top
    expected_width, expected_height = expected_size
    frame_left, frame_right = clipped_axis_window(
        left, visible_width, image_size[0], expected_width, "horizontal"
    )
    frame_top, frame_bottom = clipped_axis_window(
        top, visible_height, image_size[1], expected_height, "vertical"
    )
    scale_x = frame_size[0] / expected_width
    scale_y = frame_size[1] / expected_height
    return (
        round(frame_left * scale_x),
        round(frame_top * scale_y),
        round(frame_right * scale_x),
        round(frame_bottom * scale_y),
    )


def set_opacity(image: Image.Image, opacity: float) -> Image.Image:
    red, green, blue, alpha = image.split()
    alpha = alpha.point(lambda value: round(value * opacity))
    return Image.merge("RGBA", (red, green, blue, alpha))


def merge_frame(
    frame_path: Path,
    image_path: Path,
    bbox: tuple[int, int, int, int],
    expected_size: tuple[int, int],
) -> Path:
    with Image.open(frame_path) as source_frame:
        frame = source_frame.convert("RGBA")
    with Image.open(image_path) as source_image:
        crop_box = frame_crop_box(bbox, source_image.size, expected_size, frame.size)
        background = source_image.crop(bbox).convert("RGBA")
    frame = frame.crop(crop_box)
    background = background.resize(frame.size, resample=Image.Resampling.LANCZOS)
    background = set_opacity(background, 0.7)
    if "umap" in frame_path.name:
        frame = set_opacity(frame, 0.7)
    result = Image.alpha_composite(
        Image.new("RGBA", frame.size, (0, 0, 0, 255)), background
    )
    result = Image.alpha_composite(result, frame)
    output_path = frame_path.with_name(f"merged_{frame_path.name}")
    result.save(output_path)
    return output_path


def main() -> None:
    args = parse_args()
    bbox = frame_bbox(args.mask)
    expected_size = (
        int(
            (
                args.y_spots_number * args.length_spot
                + (args.y_spots_number - 1) * args.interval
            )
            / args.pixel_length
        ),
        int(
            (
                args.x_spots_number * args.length_spot
                + (args.x_spots_number - 1) * args.interval
            )
            / args.pixel_length
        ),
    )
    for frame_path in args.frame:
        if not frame_path.is_file():
            print(f"Skipping missing frame: {frame_path}")
            continue
        output_path = merge_frame(frame_path, args.image, bbox, expected_size)
        print(f"Merged: {frame_path} -> {output_path}")


if __name__ == "__main__":
    main()
