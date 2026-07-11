import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def saturation_model(fraction, maximum, rate):
    return maximum * (1 - np.exp(-rate * fraction))


def fit_saturation(fractions, values):
    if np.all(values == 0):
        return 0.0, 0.0, 1.0
    observed_maximum = max(values)
    parameters, _ = curve_fit(
        saturation_model,
        fractions,
        values,
        p0=(observed_maximum * 1.2, 3.0),
        bounds=(
            (observed_maximum, 1e-8),
            (observed_maximum * 100, 1000),
        ),
        maxfev=20000,
    )
    fitted = saturation_model(fractions, *parameters)
    residual_sum = np.sum((values - fitted) ** 2)
    total_sum = np.sum((values - np.mean(values)) ** 2)
    r_squared = 1 - residual_sum / total_sum if total_sum > 0 else 1.0
    return parameters[0], parameters[1], r_squared


def summarize_fraction(fraction, data_files, source):
    data = pd.concat(
        [pd.read_csv(path, usecols=["umi_count", "gene_count"]) for path in data_files],
        ignore_index=True,
    )
    return {
        "fraction": fraction,
        "source": source,
        "spots": len(data),
        "median_umi_per_spot": data["umi_count"].median(),
        "median_genes_per_spot": data["gene_count"].median(),
    }


def find_mrna_data(results_dir):
    results_dir = Path(results_dir)
    direct = results_dir / "Solo.out/GeneFull/raw/data.csv"
    if direct.is_file():
        return [direct]
    else:
        raise FileNotFoundError(f"No mRNA data.csv found in {results_dir}")


def main():
    parser = argparse.ArgumentParser(description="Plot mRNA saturation curves")
    parser.add_argument("--input", required=True, help="Saturation result directory")
    parser.add_argument("--original", required=True, help="Original mRNA results directory")
    parser.add_argument("--fractions", required=True, help="Comma-separated sampling fractions")
    parser.add_argument("--output", required=True, help="Output PNG path")
    args = parser.parse_args()

    input_dir = Path(args.input)
    rows = {}
    for value in args.fractions.split(","):
        fraction = float(value.strip())
        if np.isclose(fraction, 1.0):
            continue
        label = f"{fraction:.2f}"
        data_files = find_mrna_data(input_dir / label / "results")
        if not data_files:
            raise FileNotFoundError(f"No mRNA data.csv found for fraction {label}")

        rows[fraction] = summarize_fraction(
            fraction, data_files, source="downsampled"
        )

    original_files = find_mrna_data(args.original)
    if not original_files:
        raise FileNotFoundError(
            f"No original mRNA data.csv found below {args.original}"
        )
    rows[1.0] = summarize_fraction(1.0, original_files, source="original")

    metrics = pd.DataFrame(rows.values()).sort_values("fraction")
    if len(metrics) < 2:
        raise ValueError("At least two fractions are required to fit a saturation curve")
    metrics.to_csv(input_dir / "saturation_metrics.csv", index=False)

    fractions = metrics["fraction"].to_numpy(dtype=float)
    curve_settings = [
        ("median_umi_per_spot", "Median UMI per Spot"),
        ("median_genes_per_spot", "Median Genes per Spot"),
    ]
    fit_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for axis, (column, ylabel) in zip(axes, curve_settings):
        values = metrics[column].to_numpy(dtype=float)
        maximum, rate, r_squared = fit_saturation(fractions, values)
        fraction_95 = -np.log(0.05) / rate if rate > 0 else max(fractions)
        current_index = int(np.argmin(np.abs(fractions - 1.0)))
        current_saturation = (
            values[current_index] / maximum * 100 if maximum > 0 else 0.0
        )
        x_limit = max(max(fractions) * 1.2, fraction_95 * 1.1)
        y_limit = max(maximum, max(values)) * 1.1
        if y_limit == 0:
            y_limit = 1
        smooth_fractions = np.linspace(0, x_limit, 500)
        fit_rows.append(
            {
                "metric": column,
                "estimated_saturation": maximum,
                "rate": rate,
                "fraction_at_95_percent_saturation": fraction_95,
                "current_saturation_percent": current_saturation,
                "r_squared": r_squared,
            }
        )
        axis.scatter(fractions * 100, values, label="Observed", zorder=3)
        axis.plot(
            smooth_fractions * 100,
            saturation_model(smooth_fractions, maximum, rate),
            label=f"Fit: max={maximum:.1f}, $R^2$={r_squared:.3f}",
        )
        axis.axhline(
            maximum,
            color="red",
            linestyle="--",
            linewidth=1.2,
            label=f"Saturation={maximum:.1f}",
        )
        axis.plot(
            [],
            [],
            linestyle="none",
            label=f"Current saturation={current_saturation:.1f}%",
        )
        axis.set_xlabel("Sequencing Fraction (%)", fontsize=10)
        axis.set_ylabel(ylabel, fontsize=10)
        axis.set_xlim(0, x_limit * 100)
        axis.set_ylim(0, y_limit)
        axis.tick_params(axis="x", labelsize=10)
        axis.grid(alpha=0.3)
        axis.legend()
    pd.DataFrame(fit_rows).to_csv(input_dir / "saturation_fit.csv", index=False)
    fig.suptitle("mRNA Sequencing Saturation", fontsize=12)
    fig.tight_layout()
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
