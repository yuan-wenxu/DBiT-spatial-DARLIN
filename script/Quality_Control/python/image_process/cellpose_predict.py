from stereo.algorithm import cell_pose as cp
import tifffile
import numpy as np
import argparse

class CellposeConfig:
    def __init__(self, top_value, number_of_top_values, dmin, dmax, step, photo_size, photo_step):
        self.top_value = top_value
        self.number_of_top_values = number_of_top_values
        self.dmin = dmin
        self.dmax = dmax
        self.step = step
        self.photo_size = photo_size
        self.photo_step = photo_step

def cellpose_pre(image_path, output_path, config):
    img = tifffile.imread(image_path)
    top_values = np.partition(img[:, :, 1].flatten(), -config.number_of_top_values)[-config.number_of_top_values:]
    if np.mean(top_values) <= config.top_value:
        pass
    else:
        cp.Cellpose(
                    img_path=image_path,
                    out_path=output_path,
                    model_type='cyto2',
                    dmin=config.dmin,   # min cell diameter
                    dmax=config.dmax,   # max cell diameter
                    step=config.step,
                    gpu=False,
                    photo_size=config.photo_size,
                    photo_step=config.photo_step
                    )
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cellpose prediction')
    parser.add_argument('-i', '--input_file', type=str, help='input image')
    parser.add_argument('-o', '--output_path', type=str, help='output image')
    parser.add_argument('-t', '--top_value', type=int, default=50, help='top value of the image')
    parser.add_argument('-n', '--number_of_top_values', type=int, default=1500, help='number of top values')
    parser.add_argument('-dmin', '--dmin', type=int, default=20, help='min cell diameter')
    parser.add_argument('-dmax', '--dmax', type=int, default=50, help='max cell diameter')
    parser.add_argument('-step', '--step', type=int, default=10, help='step size')
    parser.add_argument('-photo_size', '--photo_size', type=int, default=170, help='photo size')
    parser.add_argument('-photo_step', '--photo_step', type=int, default=170, help='photo step size')
    args = parser.parse_args()
    input_file = args.input_file
    output_path = args.output_path
    config = CellposeConfig(args.top_value, args.number_of_top_values, args.dmin, args.dmax, args.step, args.photo_size, args.photo_step)
    cellpose_pre(input_file, output_path, config)