"""Coarse whole-image tissue masking for DBiT spot classification."""

from pathlib import Path

import cv2
import numpy as np


BACKGROUND_THRESHOLD = 8
LOCAL_DENSITY_KERNEL = 201
MIN_LOCAL_SIGNAL_FRACTION = 0.02
MORPHOLOGY_KERNEL = 11
MIN_COMPONENT_FRACTION = 0.0001
MIN_COMPONENT_RELATIVE_TO_LARGEST = 0.001
MIN_TISSUE_FRACTION = 0.03


def _grayscale(image):
    if image.ndim == 2:
        return image
    if image.ndim != 3:
        raise ValueError(f"Unsupported image shape: {image.shape}")
    if image.shape[2] == 1:
        return image[:, :, 0]
    return cv2.cvtColor(image[:, :, :3], cv2.COLOR_BGR2GRAY)


def _to_uint8(image):
    if image.dtype == np.uint8:
        return image
    normalized = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX)
    return normalized.astype(np.uint8)


def _remove_small_regions(mask):
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


def generate_tissue_mask(image_path, output_path):
    """Generate and save a coarse binary mask; return it as a uint8 array."""
    image = cv2.imread(str(image_path), cv2.IMREAD_UNCHANGED)
    if image is None:
        raise ValueError(f"Unable to read image: {image_path}")

    intensity = _to_uint8(_grayscale(image))
    del image
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

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(output_path), mask):
        raise RuntimeError(f"Failed to write tissue mask: {output_path}")
    return mask
