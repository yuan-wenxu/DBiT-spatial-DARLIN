import cv2
import argparse
from .stardist_predict import PlotConfig
from .tissue_mask import MIN_TISSUE_FRACTION
from tqdm import tqdm

ORIENTATION_CHOICES = ('normal', 'horizontal', 'vertical', 'rotate')


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


def normalize_orientation(value):
    orientation = value.lower()
    if orientation not in ORIENTATION_CHOICES:
        raise argparse.ArgumentTypeError(
            f'orientation must be one of {", ".join(ORIENTATION_CHOICES)}, got: {value}'
        )
    return orientation


def split_image(image_path, config, output_path, result_path, put_text, font_size,
                orientation='normal', swap_xy=False, tissue_mask=None, write_tiles=True):
    img = cv2.imread(image_path)
    if img is None:
        raise ValueError(f'Unable to read image: {image_path}')
    x_start_list = [int(x * (config.length_spot + config.interval) / config.pixel_length) for x in range(0, config.x_spots_number, 1)]
    y_start_list = [int(y * (config.length_spot + config.interval) / config.pixel_length) for y in range(0, config.y_spots_number, 1)]
    x_end_list = [x_start_list[i] + int(config.length_spot / config.pixel_length) for i in range(len(x_start_list))]
    y_end_list = [y_start_list[i] + int(config.length_spot / config.pixel_length) for i in range(len(y_start_list))]
    required_width = max(x_end_list)
    required_height = max(y_end_list)
    image_height, image_width = img.shape[:2]
    if image_width < required_width or image_height < required_height:
        raise ValueError(
            f'Image is smaller than the configured DBiT grid: '
            f'image={image_width}x{image_height}, required={required_width}x{required_height}'
        )
    if tissue_mask is not None and tissue_mask.shape[:2] != img.shape[:2]:
        raise ValueError(
            f'Tissue mask shape {tissue_mask.shape[:2]} does not match image shape {img.shape[:2]}'
        )

    output_x_count = config.y_spots_number if swap_xy else config.x_spots_number
    output_y_count = config.x_spots_number if swap_xy else config.y_spots_number
    reverse_x = orientation in ('horizontal', 'rotate')
    reverse_y = orientation in ('vertical', 'rotate')
    tissue_spots = {}

    def get_crop_indices(x, y):
        mapped_x = output_x_count - 1 - x if reverse_x else x
        mapped_y = output_y_count - 1 - y if reverse_y else y
        if swap_xy:
            return mapped_y, mapped_x
        return mapped_x, mapped_y

    description = 'Splitting image and classifying tissue' if write_tiles else 'Classifying tissue'
    for x in tqdm(range(output_x_count), desc=description):
        for y in range(output_y_count):
            crop_x, crop_y = get_crop_indices(x, y)
            x_s = x_start_list[crop_x]
            x_e = x_end_list[crop_x]
            y_s = y_start_list[crop_y]
            y_e = y_end_list[crop_y]
            start_point = (int(x_s), int(y_s))
            end_point = (int(x_e), int(y_e))
            if tissue_mask is not None:
                mask_crop = tissue_mask[start_point[1]:end_point[1], start_point[0]:end_point[0]]
                tissue_fraction = float((mask_crop > 0).sum()) / float(mask_crop.size)
                tissue_spots[(x, y)] = tissue_fraction >= MIN_TISSUE_FRACTION
            if write_tiles:
                img_small = img[start_point[1]:end_point[1], start_point[0]:end_point[0]]
                cv2.imwrite(f'{output_path}/{x}_{y}.tif', img_small)
                cv2.rectangle(img, start_point, end_point, (255, 0, 0), 2)
                if put_text:
                    cv2.putText(img, str(f'{x}_{y}'), (int(x_s), int(y_e)), cv2.FONT_HERSHEY_SIMPLEX, font_size, (0, 0, 255), 1)
    if write_tiles:
        cv2.imwrite(f'{result_path}/result.png', img)
    return tissue_spots

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Segmentation')
    parser.add_argument('-f', '--file_path', type=str, help='path to the image')
    parser.add_argument('-x', '--x_spots_number', type=int, default=50, help='Number of spots in x direction')
    parser.add_argument('-y', '--y_spots_number', type=int, default=50, help='Number of spots in y direction')
    parser.add_argument('-l', '--length_spot', type=int, default=50, help='Length of each spot in μm')
    parser.add_argument('-i', '--interval', type=int, default=50, help='Interval between two adjacent spots in μm')
    parser.add_argument('-p', '--pixel_length', type=float, default=0.294, help='Length of each pixel in μm')
    parser.add_argument('-o', '--output_path', help='path to the output image')
    parser.add_argument('-r', '--result_path', help='path to the result image')
    parser.add_argument('-t', '--put_text', type=str_to_bool, default=True, help='Whether to put text on the image (True/False, yes/no, 1/0)')
    parser.add_argument('-fs', '--font_size', type=int, default=1, help='Font size of the text')
    parser.add_argument('--orientation', type=normalize_orientation, default='normal', help='Grid origin orientation: normal, horizontal, vertical, or rotate')
    parser.add_argument('--swap_xy', action='store_true', help='Swap x and y grid axes after applying orientation')
    args = parser.parse_args()

    image_path = args.file_path
    x_spots_number = args.x_spots_number
    y_spots_number = args.y_spots_number
    length_spot = args.length_spot
    interval = args.interval
    pixel_length = args.pixel_length
    output_path = args.output_path
    result_path = args.result_path
    put_text = args.put_text
    font_size = args.font_size
    orientation = args.orientation
    swap_xy = args.swap_xy

    config = PlotConfig(x_spots_number, y_spots_number, length_spot, interval, pixel_length)
    split_image(image_path, config, output_path, result_path, put_text=put_text, font_size=font_size, orientation=orientation, swap_xy=swap_xy)
