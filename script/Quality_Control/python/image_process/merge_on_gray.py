"""
Vertically flip *_filtered.png images and composite them onto gray.png.

Usage:
    python merge_on_gray.py --gray <gray.png> --search-dir <dir> [--recursive]

For each *_filtered.png found under --search-dir the script will:
  1. Flip the image vertically.
  2. Resize gray.png to match the flipped image.
  3. Composite the RGBA overlay onto the gray background.
  4. Save as merged_<original_filename> in the same directory.
"""
import argparse
import glob
import os
from pathlib import Path

from PIL import Image, ImageOps
Image.MAX_IMAGE_PIXELS = None 

def set_opacity(img: Image.Image, opacity: float) -> Image.Image:
    """Return a copy of img with its alpha channel scaled by opacity (0-1)."""
    r, g, b, a = img.split()
    a = a.point(lambda x: round(x * opacity))
    return Image.merge("RGBA", (r, g, b, a))


def merge_on_gray(frame_path: str, gray_path: str) -> None:
    frame = Image.open(frame_path).convert("RGBA")
    if "umap" in frame_path:
        frame_flipped = set_opacity(ImageOps.flip(frame), 0.7)
    else:
        frame_flipped = ImageOps.flip(frame)

    gray = Image.open(gray_path).convert("RGBA")
    gray_resized = set_opacity(gray.resize(frame_flipped.size, resample=Image.LANCZOS), 0.7)

    # Black opaque background ensures the final image is fully opaque
    black_bg = Image.new("RGBA", frame_flipped.size, (0, 0, 0, 255))
    result = Image.alpha_composite(black_bg, gray_resized)
    result = Image.alpha_composite(result, frame_flipped)

    out_dir = os.path.dirname(frame_path)
    out_name = "merged_" + os.path.basename(frame_path)
    result.save(os.path.join(out_dir, out_name))


def main() -> None:
    parser = argparse.ArgumentParser(description="Flip and merge frame plots onto gray background.")
    parser.add_argument("--gray", required=True, help="Path to gray.png background image")
    parser.add_argument("--search-dir", required=True, help="Directory to search for *_filtered.png")
    parser.add_argument("--recursive", action="store_true", help="Search subdirectories recursively")
    args = parser.parse_args()

    if not os.path.isfile(args.gray):
        raise FileNotFoundError(f"gray image not found: {args.gray}")

    pattern = "*_filtered.png"
    if args.recursive:
        matches = list(Path(args.search_dir).rglob(pattern))
    else:
        matches = list(Path(args.search_dir).glob(pattern))

    if not matches:
        print(f"No {pattern} files found under {args.search_dir}")
        return

    for p in sorted(matches):
        merge_on_gray(str(p), args.gray)
        print(f"Merged: {p} -> merged_{p.name}")


if __name__ == "__main__":
    main()
