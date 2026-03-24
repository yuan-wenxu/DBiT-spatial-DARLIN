from image_process import split_image, cell_num_area, cellpose_pre
from image_process import CellposeConfig
from plot import PlotConfig
import os
import pandas as pd
import argparse

def str_to_bool(value):
    """Convert a string to boolean."""
    if isinstance(value, bool):
        return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError(f'Boolean value expected, got: {value}')

def main(image_path, output_path, result_path, label_path, mask_path, config, cellpose_config, put_text, font_size):
    result = pd.DataFrame(columns=['x', 'y', 'num_cells', 'area', 'status'])
    put_text = str_to_bool(put_text)
    split_image(image_path, config, output_path, result_path, put_text, font_size)
    for i in os.listdir(output_path):
        mask_file = os.path.join(mask_path, i)
        cellpose_pre(os.path.join(output_path, i), mask_file, cellpose_config)
        if os.path.exists(mask_file):
            label_file = os.path.join(label_path, i)
            df = cell_num_area(os.path.join(output_path, i), mask_file, label_file)
            result = result.append(df, ignore_index=True)
        else:
            x = i.split('.')[0].split('_')[0]
            y = i.split('.')[0].split('_')[1]
            result = result.append({'x': x, 'y': y, 'num_cells': 0, 'area': [],'status': 'skiped'}, ignore_index=True)
    result.to_csv(os.path.join(result_path, 'cell_num_area.csv'), index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process image and get cell number and area')
    parser.add_argument('-ip', '--image_path', type=str, help='path of image file')
    parser.add_argument('-r', '--result_path', type=str, help='path of result file')
    parser.add_argument('-x', '--x_spots_number', type=int, default=50, help='Number of spots in x direction')
    parser.add_argument('-y', '--y_spots_number', type=int, default=50, help='Number of spots in y direction')
    parser.add_argument('-l', '--length_spot', type=int, default=50, help='Length of each spot in μm')
    parser.add_argument('-i', '--interval', type=int, default=50, help='Interval between two adjacent spots in μm')
    parser.add_argument('-p', '--pixel_length', type=float, default=0.294, help='Length of each pixel in μm')
    parser.add_argument('-t', '--put_text', type=str_to_bool, default=True, help='Whether to put text on the image (True/False, yes/no, 1/0)')
    parser.add_argument('-fs', '--font_size', type=int, default=1, help='Font size of the text')
    parser.add_argument('-top_value', '--top_value', type=int, default=50, help='top value of the image')
    parser.add_argument('-number_of_top_values', '--number_of_top_values', type=int, default=1500, help='number of top values')
    parser.add_argument('-dmin', '--dmin', type=int, default=20, help='min cell diameter')
    parser.add_argument('-dmax', '--dmax', type=int, default=50, help='max cell diameter')
    parser.add_argument('-step', '--step', type=int, default=10, help='step size')
    parser.add_argument('-photo_size', '--photo_size', type=int, default=170, help='photo size')
    parser.add_argument('-photo_step', '--photo_step', type=int, default=170, help='photo step size')
    args = parser.parse_args()

    image_path = args.image_path
    result_path = args.result_path
    output_path = result_path + '/' + 'split'
    label_path = result_path + '/' + 'label'
    mask_path = result_path + '/' + 'mask'
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(label_path, exist_ok=True)
    os.makedirs(mask_path, exist_ok=True)
    
    x_spots_number = args.x_spots_number
    y_spots_number = args.y_spots_number
    length_spot = args.length_spot
    interval = args.interval
    pixel_length = args.pixel_length
    config = PlotConfig(x_spots_number, y_spots_number, length_spot, interval, pixel_length)
    top_value = args.top_value
    number_of_top_values = args.number_of_top_values
    dmin = args.dmin
    dmax = args.dmax
    step = args.step
    photo_size = args.photo_size
    photo_step = args.photo_step
    cellpose_config = CellposeConfig(top_value, number_of_top_values, dmin, dmax, step, photo_size, photo_step)
    put_text = args.put_text
    font_size = args.font_size
    main(image_path, output_path, result_path, label_path, mask_path, config, cellpose_config, put_text, font_size)