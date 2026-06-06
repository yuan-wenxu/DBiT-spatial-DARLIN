from image_process import stardist_pre, StardistConfig, split_image, PlotConfig
import os
import pandas as pd
import argparse
import numpy as np


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


ORIENTATION_CHOICES = ('normal', 'horizontal', 'vertical', 'rotate')


def normalize_orientation(value):
    orientation = value.lower()
    if orientation not in ORIENTATION_CHOICES:
        raise argparse.ArgumentTypeError(
            f'orientation must be one of {", ".join(ORIENTATION_CHOICES)}, got: {value}'
        )
    return orientation


def main(image_path, split_path, mask_path, label_path, plot_config, stardist_config, put_text, font_size, orientation):
    os.makedirs(split_path, exist_ok=True)
    os.makedirs(mask_path, exist_ok=True)
    os.makedirs(label_path, exist_ok=True)
    result_path = os.path.dirname(split_path)
    orientation_file = os.path.join(split_path, '.orientation')
    expected_split_count = plot_config.x_spots_number * plot_config.y_spots_number
    split_files = [
        i for i in os.listdir(split_path)
        if i.lower().endswith(('.tif', '.tiff', '.png', '.jpg'))
    ]
    previous_orientation = ''
    if os.path.exists(orientation_file):
        with open(orientation_file, 'r', encoding='utf-8') as f:
            previous_orientation = f.read().strip()

    if len(split_files) == expected_split_count and previous_orientation == orientation:
        print(f'✓ Split images already exist in {split_path} with orientation={orientation}, skipping splitting step.')
    else:
        split_image(image_path, plot_config, split_path, result_path, put_text, font_size, orientation)
        with open(orientation_file, 'w', encoding='utf-8') as f:
            f.write(orientation)

    result = pd.DataFrame(columns=['x', 'y', 'num_cells', 'area', 'status'])

    for i in os.listdir(split_path):
        if not i.lower().endswith(('.tif', '.tiff', '.png', '.jpg')):
            continue
        
        split_file = os.path.join(split_path, i)
        mask_file = os.path.join(mask_path, i)
        label_file = os.path.join(label_path, i)
        
        try:
            labels = stardist_pre(split_file, mask_file, label_file, stardist_config)
            
            if labels is not None:
                # Extract x, y coordinates from filename
                x = i.split('.')[0].split('_')[0]
                y = i.split('.')[0].split('_')[1]
                num_cells = int(np.max(labels))
                
                # Get cell areas
                areas = []
                for cell_id in range(1, num_cells + 1):
                    # Cast to native Python int to avoid numpy int64 representation in CSV
                    area = int(np.sum(labels == cell_id))
                    areas.append(area)
                
                result = pd.concat(
                    [result, pd.DataFrame({
                        'x': [x], 
                        'y': [y], 
                        'num_cells': [num_cells], 
                        'area': [areas], 
                        'status': ['predicted']
                    })],
                    ignore_index=True
                )
                print(f'✓ Processed: {i} - {num_cells} cells')
            else:
                # Low quality image - skip
                x = i.split('.')[0].split('_')[0]
                y = i.split('.')[0].split('_')[1]
                result = pd.concat(
                    [result, pd.DataFrame({
                        'x': [x], 
                        'y': [y], 
                        'num_cells': [0], 
                        'area': [[]], 
                        'status': ['skipped']
                    })],
                    ignore_index=True
                )
                print(f'✗ Skipped: {i} - no cells')
        except Exception as e:
            x = i.split('.')[0].split('_')[0]
            y = i.split('.')[0].split('_')[1]
            result = pd.concat(
                [result, pd.DataFrame({
                    'x': [x], 
                    'y': [y], 
                    'num_cells': [0], 
                    'area': [[]], 
                    'status': [f'error: {str(e)}']
                })],
                ignore_index=True
            )
            print(f'✗ Error processing {i}: {str(e)}')

    result.to_csv(os.path.join(result_path, 'cell_num_area.csv'), index=False)
    print(f'\n=== Processing Complete ===')
    print(f'Total processed: {len(result)}')
    print(f'Successfully predicted: {len(result[result["status"] == "predicted"])}')
    print(f'Skipped/Failed: {len(result[result["status"] != "predicted"])}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='StarDist cell segmentation')
    parser.add_argument('-ip', '--image_path', type=str, help='path of split image folder')
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
    parser.add_argument('-m', '--model_name', type=str, default='2D_versatile_fluo', help='pretrained model name')
    parser.add_argument('-pt', '--prob_thresh', type=float, default=0.5, help='detection probability threshold')
    parser.add_argument('-nt', '--nms_thresh', type=float, default=0.6, help='NMS IoU threshold')
    parser.add_argument('--orientation', type=normalize_orientation, default='normal', help='Grid origin orientation: normal, horizontal, vertical, or rotate')
    args = parser.parse_args()

    result_path = args.result_path
    split_path = result_path + '/' + 'split'
    mask_path = result_path + '/' + 'mask'
    label_path = result_path + '/' + 'label'

    plot_config = PlotConfig(
        args.x_spots_number, args.y_spots_number, args.length_spot, args.interval, args.pixel_length
    )

    stardist_config = StardistConfig(
        args.top_value, args.number_of_top_values, args.prob_thresh, args.nms_thresh, args.model_name
    )
    main(args.image_path, split_path, mask_path, label_path, plot_config, stardist_config, args.put_text, args.font_size, args.orientation)
