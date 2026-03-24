import tifffile
import numpy as np
from skimage import measure
from skimage.measure import regionprops
from skimage.segmentation import mark_boundaries
from skimage.io import imsave
import os
import pandas as pd
import argparse

def cell(image_path, mask_path, label_path):
    mask = tifffile.imread(mask_path)
    binary_mask = (mask > 0).astype(np.uint8)
    label_image = measure.label(binary_mask, connectivity=1)
    num_cells = len(np.unique(label_image)) - 1
    regions = regionprops(label_image)
    areas = [region.area for region in regions]
    image_name = os.path.basename(image_path)
    x = image_name.split('.')[0].split('_')[0]
    y = image_name.split('.')[0].split('_')[1]
    img = tifffile.imread(image_path)
    overlay = mark_boundaries(img, label_image, color=(1, 0, 0))
    overlay = (overlay * 255).astype(np.uint8)
    imsave(label_path, overlay)
    return pd.DataFrame({'x': x, 'y': y, 'num_cells': num_cells, 'area': [areas], 'status': 'predicted'})

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Cell Counting')
    parser.add_argument('-i', '--image_path', type=str, help='path to image file')
    parser.add_argument('-m', '--mask_path', type=str, help='path to mask file')
    parser.add_argument('-l', '--label_path', type=str, help='path to label file')
    args = parser.parse_args()
    df = cell(args.image_path, args.mask_path, args.label_path)