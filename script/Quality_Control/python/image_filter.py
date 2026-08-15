#!/usr/bin/env python3
"""Filter segmented cells by area and summarize retained cells per DBiT spot."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def filter_segmented_cells(result_path: Path, cutoff: float) -> pd.DataFrame:
    cell_file = result_path / "cell_num_area.csv"
    area_plot = result_path / "area.png"
    filtered_file = result_path / "filtered_results.csv"

    cells = pd.read_csv(cell_file)
    cells["area"] = cells["area"].apply(ast.literal_eval)
    nonempty_areas = [np.asarray(values) for values in cells["area"] if values]
    all_areas = (
        np.concatenate(nonempty_areas) if nonempty_areas else np.array([], dtype=float)
    )

    print(f"Cell number before filtering: {len(all_areas)}")
    print(f'Spots number before filtering: {(cells["num_cells"] != 0).sum()}')

    figure, axis = plt.subplots(figsize=(5, 4))
    axis.hist(all_areas, bins=100)
    axis.axvline(
        cutoff, linestyle="--", linewidth=1, color="r", label=f"cutoff = {cutoff}"
    )
    axis.set_xlabel("Area")
    axis.set_ylabel("Count")
    axis.set_title("Cell area distribution")
    axis.legend()
    figure.tight_layout()
    figure.savefig(area_plot, dpi=300, bbox_inches="tight")
    plt.close(figure)

    cells["count"] = [
        sum(area >= cutoff for area in values) for values in cells["area"]
    ]
    print(f'Cell number after filtering: {cells["count"].sum()}')
    print(f'Spots number after filtering: {(cells["count"] != 0).sum()}')
    cells.to_csv(filtered_file, index=False)
    return cells


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter StarDist cells by area and count retained cells per spot."
    )
    parser.add_argument(
        "--file_path",
        required=True,
        type=Path,
        help="Directory containing cell_num_area.csv",
    )
    parser.add_argument(
        "--cutoff", required=True, type=float, help="Minimum retained area"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cells = filter_segmented_cells(args.file_path, args.cutoff)


if __name__ == "__main__":
    main()
