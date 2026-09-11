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


def set_opacity(image: Image.Image, opacity: float) -> Image.Image:
    red, green, blue, alpha = image.split()
    alpha = alpha.point(lambda value: round(value * opacity))
    return Image.merge("RGBA", (red, green, blue, alpha))


def merge_frame(
    frame_path: Path,
    image_path: Path,
    bbox: tuple[int, int, int, int],
) -> Path:
    with Image.open(frame_path) as source_frame:
        frame = source_frame.convert("RGBA")
    with Image.open(image_path) as source_image:
        background = source_image.crop(bbox).convert("RGBA")
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
    for frame_path in args.frame:
        if not frame_path.is_file():
            print(f"Skipping missing frame: {frame_path}")
            continue
        output_path = merge_frame(frame_path, args.image, bbox)
        print(f"Merged: {frame_path} -> {output_path}")


if __name__ == "__main__":
    main()
