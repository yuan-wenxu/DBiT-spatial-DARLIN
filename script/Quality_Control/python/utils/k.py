import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick

def k_category(k):
    if k <= 1:
        return "≤1"
    elif k <= 5:
        return "≤5"
    elif k <= 10:
        return "≤10"
    else:
        return ">10"
    
def k_disperse(df, output):
    df_summary = (df.groupby("SR").agg(n_reads=("reads", "sum"), n_UR=("UR", "nunique")).reset_index())
    df_summary["k"] = df_summary["n_reads"] / df_summary["n_UR"]
    df_summary["k_category"] = df_summary["k"].apply(k_category)
    k_cat_order = ["≤1", "≤5", "≤10", ">10"]
    k_cat_colors = {
        "≤1": "#4575b4",
        "≤5": "#91bfdb",
        "≤10": "#fee090",
        ">10": "#d73027"
        }
    plt.figure(figsize=(5, 3))
    for cat in k_cat_order:
        data = df_summary[df_summary["k_category"] == cat]
        plt.scatter(
            data["n_reads"],
            data["n_UR"],
            s=2,
            alpha=0.4,
            color=k_cat_colors[cat],
            label=cat
        )
    if np.max(data["n_reads"]) >= 10:
        plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Reads")
    plt.ylabel("UMIs")
    x_min = df_summary["n_UR"].min()
    x_max = df_summary["n_UR"].max()
    x_vals = np.array([x_min, x_max])
    plt.plot(x_vals, x_vals, linestyle="--", color="red", linewidth=1, label="slope=1")
    legend_handles = [mpatches.Patch(color=k_cat_colors[cat], label=f"k {cat}") for cat in k_cat_order]
    plt.legend(handles=legend_handles, title="k = Reads/UMIs", loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f"{output}/k.png")

    upper = 10
    k_vals = df_summary["k"].clip(upper=upper)
    plt.figure(figsize=(4, 3))
    plt.hist(k_vals, bins=range(1, upper + 2), color="#4575b4", alpha=0.7, edgecolor="black")
    plt.xlabel(f"Reads / UMIs (k, capped at {upper})")
    plt.ylabel("Spots number")
    plt.title("Histogram of k")
    plt.xticks(list(range(1, upper + 1)) + [upper], labels=[str(i) for i in range(1, upper + 1)] + [f"{upper}+"])
    plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(f"{output}/k_hist.png")

    k_cutoffs = range(1, 16)
    n_SR_above_k = [np.sum(df_summary["k"] >= k) for k in k_cutoffs]
    plt.figure(figsize=(4, 3))
    plt.plot(k_cutoffs, n_SR_above_k, marker="o", color="#4575b4", alpha=0.8)
    plt.xlabel("k cutoff (Reads / UMIs)")
    plt.ylabel("Number of SRs ≥ k cutoff")
    plt.gca().yaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # plt.yscale("log")
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(f"{output}/k_cutoff_scatter.png")

    return df_summary