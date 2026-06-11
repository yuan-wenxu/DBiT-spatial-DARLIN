import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


def setup_plot_dir(args):
    args.output_path.mkdir(parents=True, exist_ok=True)
    return args.output_path


def finish_plot(fig, out_file):
    for ax in fig.axes:
        ax.xaxis.label.set_size(10)
        ax.yaxis.label.set_size(10)
        ax.tick_params(axis="both", labelsize=10)
        ax.title.set_size(12)
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote plot: {out_file}")


def plot_lineage_length_hist(lb_len, out_file, min_len=None):
    if len(lb_len) == 0:
        return
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(lb_len, bins=range(1, 300, 1), edgecolor="black")
    if min_len is not None:
        ax.axvline(min_len, color="red", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Sequence Length", fontsize=10)
    ax.set_ylabel("Number of reads", fontsize=10)
    ax.set_title("Distribution of DARLIN Array Sequence\nLengths By Reads", fontsize=12)
    finish_plot(fig, out_file)


def get_cutoff_values(df):
    if df.empty:
        return []
    reads = df["reads"]
    max_reads = int(reads.max())
    if max_reads < 1:
        return []
    values = list(range(1, min(10, max_reads) + 1))
    if max_reads >= 11:
        values.extend(range(11, min(50, max_reads) + 1, 3))
    coarse_stop = max_reads // 2
    if max_reads >= 61 and coarse_stop >= 61:
        values.extend(range(61, coarse_stop + 1, 10))
    return sorted(set(values))


def plot_reads_cutoff_qc(df, reads_cutoff, out_file):
    cutoff_values = get_cutoff_values(df)
    if not cutoff_values:
        return
    total_reads = df["reads"].sum()
    num_molecules = [df.loc[df["reads"] >= c, "UB"].nunique() for c in cutoff_values]
    fraction_reads_retained = [df.loc[df["reads"] >= c, "reads"].sum() / total_reads for c in cutoff_values]

    fig, axes = plt.subplots(2, 1, figsize=(5, 5))
    axes[0].plot(cutoff_values, num_molecules, marker="o", markersize=2, linewidth=1)
    axes[0].axvline(reads_cutoff, color="red", linestyle="--", linewidth=0.6)
    axes[0].set_xlabel("Reads Cutoff", fontsize=10)
    axes[0].set_ylabel("Number of Molecules", fontsize=10)
    axes[0].set_yscale("log")
    axes[0].set_xscale("log")
    axes[0].grid(alpha=0.3)

    axes[1].plot(cutoff_values, fraction_reads_retained, marker="o", markersize=2, linewidth=1)
    axes[1].axvline(reads_cutoff, color="red", linestyle="--", linewidth=0.6)
    axes[1].set_xlabel("Reads Cutoff", fontsize=10)
    axes[1].set_ylabel("Frac. of Reads\nRetained", fontsize=10)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xscale("log")
    axes[1].grid(alpha=0.3)
    finish_plot(fig, out_file)


def plot_reads_fraction_qc(df, threshold, out_file):
    if df.empty:
        return
    fig, axes = plt.subplots(2, 1, figsize=(4, 5))
    axes[0].hist(df["reads_fraction"], bins=50, edgecolor="white", linewidth=0.3)
    axes[0].axvline(threshold, color="red", linestyle="--", linewidth=0.8)
    axes[0].set_xlabel("Reads Fraction", fontsize=10)
    axes[0].set_ylabel("Number of (SR, UR, LR)", fontsize=10)

    axes[1].scatter(df["reads_fraction"], df["reads"], s=0.2, alpha=0.15)
    axes[1].axvline(threshold, color="red", linestyle="--", linewidth=0.8)
    axes[1].set_xlabel("Reads Fraction", fontsize=10)
    axes[1].set_ylabel("Reads", fontsize=10)
    axes[1].set_yscale("log")
    finish_plot(fig, out_file)


def plot_sr_reads_umis(sr_summary, out_file):
    if sr_summary.empty:
        return
    fig, ax = plt.subplots(figsize=(4, 3))
    categories = [
        ("<=1", sr_summary["k"] <= 1, "#4575b4"),
        ("<=5", (sr_summary["k"] > 1) & (sr_summary["k"] <= 5), "#91bfdb"),
        ("<=10", (sr_summary["k"] > 5) & (sr_summary["k"] <= 10), "#fee090"),
        (">10", sr_summary["k"] > 10, "#d73027"),
    ]
    for _label, mask, color in categories:
        data = sr_summary.loc[mask]
        ax.scatter(data["n_reads"], data["n_UR"], s=2, alpha=0.4, color=color)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Reads")
    ax.set_ylabel("UMIs")

    x_min = max(1, sr_summary["n_UR"].min())
    x_max = sr_summary["n_UR"].max()
    if x_min < x_max:
        ax.plot([x_min, x_max], [x_min, x_max], linestyle="--", color="red", linewidth=1, label="slope=1")

    legend_handles = [mpatches.Patch(color=color, label=f"k {label}") for label, _mask, color in categories]
    ax.legend(handles=legend_handles, title="k = Reads/UMIs", loc="center left", bbox_to_anchor=(1, 0.5), fontsize=8, title_fontsize=9)
    finish_plot(fig, out_file)


def plot_lr_per_sr(df_final, out_file):
    if df_final.empty or "n_LR" not in df_final.columns:
        return
    df_plot = df_final[["SR", "n_LR"]].drop_duplicates()
    max_lr = max(8, int(df_plot["n_LR"].max()) + 2)
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.hist(df_plot["n_LR"], bins=range(1, max_lr), edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Number of LRs per SR")
    ax.set_ylabel("Number of SRs")
    ax.set_yscale("log")
    finish_plot(fig, out_file)
