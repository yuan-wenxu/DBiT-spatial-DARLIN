import io
import sys
import numpy as np
import tifffile
import argparse
from stardist.models import StarDist2D
from skimage.segmentation import mark_boundaries


class StardistConfig:
    def __init__(self, top_value, number_of_top_values, prob_thresh, nms_thresh, model_name):
        self.top_value = top_value
        self.number_of_top_values = number_of_top_values
        self.prob_thresh = prob_thresh
        self.nms_thresh = nms_thresh
        self.model_name = model_name


class PlotConfig:
    def __init__(self, x_spots_number, y_spots_number, length_spot, interval, pixel_length):
        self.x_spots_number = x_spots_number
        self.y_spots_number = y_spots_number
        self.length_spot = length_spot
        self.interval = interval
        self.pixel_length = pixel_length


_model_cache = {}


def _get_model(model_name):
    if model_name not in _model_cache:
        _stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            _model_cache[model_name] = StarDist2D.from_pretrained(model_name)
        finally:
            sys.stdout = _stdout
    return _model_cache[model_name]


def stardist_pre(image_path, output_path, label_patn, config):
    img_original = tifffile.imread(image_path)
    img = img_original.copy()

    # Intensity filter: check green channel brightness
    if img.ndim == 3 and img.shape[-1] >= 2:
        top_values = np.partition(
            img[:, :, 1].flatten(), -config.number_of_top_values
        )[-config.number_of_top_values:]
    else:
        top_values = np.partition(
            img.flatten(), -config.number_of_top_values
        )[-config.number_of_top_values:]
    
    if np.mean(top_values) <= config.top_value:
        return None

    if img.ndim == 3 and img.shape[-1] == 3:
        img_normalized = img[:, :, 1]
    elif img.ndim == 3 and img.shape[-1] == 1:
        img_normalized = img[:, :, 0]
    else:
        img_normalized = img.copy()

    img_normalized = img_normalized.astype(np.float32)
    if img_normalized.max() > 0:
        img_normalized /= img_normalized.max()

    model = _get_model(config.model_name)
    labels, _ = model.predict_instances(
        img_normalized,
        prob_thresh=config.prob_thresh,
        nms_thresh=config.nms_thresh
    )
    
    # Save visualization with black background and white mask
    binary_mask = (labels > 0).astype(np.uint8) * 255  # 黑色背景(0)，白色mask(255)
    tifffile.imwrite(output_path, binary_mask)
    
    # Mark cell boundaries on original image
    # Normalize original image for visualization
    if img_original.ndim == 3 and img_original.shape[-1] == 3:
        img_vis = img_original.copy().astype(np.float32)
        img_vis = img_vis / img_vis.max() if img_vis.max() > 0 else img_vis
    elif img_original.ndim == 3 and img_original.shape[-1] == 1:
        green = img_original[:, :, 0].astype(np.float32)
        zeros = np.zeros_like(green, dtype=np.float32)
        img_vis = np.stack([zeros, green, zeros], axis=2)
        img_vis = img_vis / img_vis.max() if img_vis.max() > 0 else img_vis
    else:
        img_vis = np.repeat(img_original[:, :, np.newaxis], 3, axis=2).astype(np.float32)
        img_vis = img_vis / img_vis.max() if img_vis.max() > 0 else img_vis
    
    # Mark boundaries with red color
    overlay = mark_boundaries(img_vis, labels, color=(1, 0, 0), mode='thick')
    overlay_uint8 = (overlay * 255).astype(np.uint8)
    tifffile.imwrite(label_patn, overlay_uint8)
    
    return labels


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='StarDist prediction')
    parser.add_argument('-i', '--input_file', type=str, help='input image')
    parser.add_argument('-o', '--output_path', type=str, help='output mask')
    parser.add_argument('-l', '--label_path', type=str, help='output label image')
    parser.add_argument('-p', '--prob_thresh', type=float, default=0.5, help='detection probability threshold')
    parser.add_argument('-n', '--nms_thresh', type=float, default=0.6, help='NMS IoU threshold')
    parser.add_argument('-m', '--model_name', type=str, default='2D_versatile_fluo', help='pretrained model name')
    parser.add_argument('-t', '--top_value', type=int, default=50, help='top brightness threshold')
    parser.add_argument('-n_top', '--number_of_top_values', type=int, default=1500, help='number of top values')
    args = parser.parse_args()

    config = StardistConfig(
        args.top_value, args.number_of_top_values, args.prob_thresh, args.nms_thresh, args.model_name
    )
    stardist_pre(args.input_file, args.output_path, args.label_path, config)