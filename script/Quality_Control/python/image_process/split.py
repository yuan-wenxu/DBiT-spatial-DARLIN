import cv2
import argparse
from .stardist_predict import PlotConfig
from tqdm import tqdm

def split_image(image_path, config, output_path, result_path, put_text, font_size):
    img = cv2.imread(image_path)
    x_start_list = [int(x * (config.length_spot + config.interval) / config.pixel_length) for x in range(0, config.x_spots_number, 1)]
    y_start_list = [int(y * (config.length_spot + config.interval) / config.pixel_length) for y in range(0, config.y_spots_number, 1)]
    x_end_list = [x_start_list[i] + int(config.length_spot / config.pixel_length) for i in range(len(x_start_list))]
    y_end_list = [y_start_list[i] + int(config.length_spot / config.pixel_length) for i in range(len(y_start_list))]
    x=0
    for x_s, x_e in tqdm(zip(x_start_list, x_end_list), desc='Splitting image'):
        y=49
        for y_s, y_e in zip(y_start_list, y_end_list):
            start_point = (int(x_s), int(y_s))
            end_point = (int(x_e), int(y_e))
            img_small = img[start_point[1]:end_point[1], start_point[0]:end_point[0]]
            cv2.imwrite(f'{output_path}/{x}_{y}.tif', img_small)
            cv2.rectangle(img, start_point, end_point, (255, 0, 0), 2)
            if put_text:
                cv2.putText(img, str(f'{x}_{y}'), (int(x_s), int(y_e)), cv2.FONT_HERSHEY_SIMPLEX, font_size, (0, 0, 255), 1)
            y-=1
        x+=1
    cv2.imwrite(f'{result_path}/result.png', img)

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
    args = parser.parse_args()

    image_path = args.file_path
    x_spots_number = args.x_spots_number
    y_spots_number = args.y_spots_number
    length_spot = args.length_spot
    interval = args.interval
    pixel_length = args.pixel_length
    output_path = args.output_path
    result_path = args.result_path

    config = PlotConfig(x_spots_number, y_spots_number, length_spot, interval, pixel_length)
    split_image(image_path, config, output_path, result_path)