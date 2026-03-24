import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import argparse

class ScatterConfig:
    def __init__(self, xlabel, ylabel, title, equal_axis, grid, saturation_x, saturation_y):
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        self.equal_axis = equal_axis
        self.grid = grid
        self.saturation_x = saturation_x
        self.saturation_y = saturation_y

def auto_rotate_xticks(p, threshold=10):
    labels = [label.get_text() for label in p.get_xticklabels()]
    if not labels:
        return
    max_len = max(len(label) for label in labels)
    if max_len > threshold or len(labels) > 10:
        plt.setp(p.get_xticklabels(), rotation=45, ha='right')
    else:
        plt.setp(p.get_xticklabels(), rotation=0, ha='center')

def plot_scatter(x, y, output, config):
    corr_p, p_val_p = pearsonr(x, y)
    corr_s, p_val_s = spearmanr(x, y)
    p = f'Pearson Correlation: {corr_p:.3f} (p={p_val_p:.3e})'
    s = f'Spearman Correlation: {corr_s:.3f} (p={p_val_s:.3e})'

    plt.figure(figsize=(5, 4))
    if config.saturation_x:
        x = x / config.saturation_x
    if config.saturation_y:
        y = y / config.saturation_y
    plt.scatter(x, y, s=10, alpha=0.1)
    auto_rotate_xticks(plt.gca())
    if config.equal_axis:
        _, xmax = plt.xlim()
        _, ymax = plt.ylim()
        max_val = max(xmax, ymax)
        plt.xlim(0, max_val)
        plt.ylim(0, max_val)
        plt.plot([0, max_val], [0, max_val], linestyle='--', color='red')
    plt.text(0.05, 0.95, p, transform=plt.gca().transAxes, fontsize=8, verticalalignment='top')
    plt.text(0.05, 0.9, s, transform=plt.gca().transAxes, fontsize=8, verticalalignment='top')
    if config.grid:
        plt.grid(True)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(config.xlabel, fontsize=14)
    plt.ylabel(config.ylabel, fontsize=14)
    plt.title(config.title, fontsize=20)
    plt.tight_layout()
    plt.savefig(f'{output}/{config.title}_scatter.png', dpi=600)

def plot_scatter_small(x, y, output, config, scope_x=40000, scope_y=40000):
    plt.figure(figsize=(3, 3))
    if config.saturation_x:
        x = x / config.saturation_x
    if config.saturation_y:
        y = y / config.saturation_y
    if scope_x:
        plt.xlim(0, scope_x)
    if scope_y:
        plt.ylim(0, scope_y)
    plt.scatter(x, y, s=10, alpha=0.1)
    auto_rotate_xticks(plt.gca(), threshold=3)
    if config.equal_axis:
        _, xmax = plt.xlim()
        _, ymax = plt.ylim()
        max_val = max(xmax, ymax)
        plt.xlim(0, max_val)
        plt.ylim(0, max_val)
        plt.plot([0, max_val], [0, max_val], linestyle='--', color='red')
    if config.grid:
        plt.grid(True)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel(config.xlabel, fontsize=10)
    plt.ylabel(config.ylabel, fontsize=10)
    plt.title(config.title, fontsize=10)
    plt.tight_layout()
    plt.savefig(f'{output}/{config.title}_scatter_small.png', dpi=600)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scatter plot')
    parser.add_argument('-x', '--x', help='x values')
    parser.add_argument('-y', '--y', help='y values')
    parser.add_argument('-o', '--output', help='output directory')
    parser.add_argument('-xl', '--xlabel', help='x label')
    parser.add_argument('-yl', '--ylabel', help='y label')
    parser.add_argument('-t', '--title', help='title')
    parser.add_argument('-ea', '--equal_axis', type=bool, help='equal axis')
    parser.add_argument('-g', '--grid', type=bool, help='grid')
    parser.add_argument('-sx', '--saturation_x', type=float, help='saturation')
    parser.add_argument('-sy', '--saturation_y', type=float, help='saturation')
    parser.add_argument('-scope_x', '--scope_x', type=float, help='scope')
    parser.add_argument('-scope_y', '--scope_y', type=float, help='scope')
    args = parser.parse_args()
    config = ScatterConfig(args.xlabel, args.ylabel, args.title, args.equal_axis, args.grid, args.saturation_x, args.saturation_y)
    plot_scatter(args.x, args.y, args.output, config)
    plot_scatter_small(args.x, args.y, args.output, config, args.scope_x, args.scope_y)