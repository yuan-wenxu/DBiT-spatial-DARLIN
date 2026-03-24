import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

def group_fraction(df, output, config):
    df['group_reads'] = (df.groupby(["SR", "UR"])["reads"].transform('sum'))
    df["reads_fraction"] = df["reads"] / df["group_reads"]
    df_major = df[df["reads_fraction"] >= config.major_fraction_threshold_molecule].copy()

    print('\n')
    print("After per-(SR, UR) major selection:")
    print('Number of reads:', df_major["reads"].sum())
    print('Number of molecules:', len(df_major))
    print("SR unique:", df_major["SR"].nunique())
    print("UR unique:", df_major["UR"].nunique())
    print("LR unique:", df_major["LR"].nunique())
    print('Reads with amplification error:', df["reads"].sum() - df_major["reads"].sum())

    fig, axes = plt.subplots(2, 1, figsize=(3, 4))
    # Histogram of reads_fraction
    axes[0].hist(df['reads_fraction'], bins=50, edgecolor='white')
    axes[0].axvline(config.major_fraction_threshold_molecule, color='red', linestyle='--', linewidth=.6)
    axes[0].set_ylabel('Number of (SR, UR)')
    axes[0].set_title("Lineage barcode filtration")
    axes[0].yaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
    axes[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # Scatter of (reads_fraction, reads)
    axes[1].scatter(df['reads_fraction'], df['reads'], s=1, alpha=0.01)
    axes[1].axvline(config.major_fraction_threshold_molecule, color='red', linestyle='--', linewidth=0.6)
    axes[1].set_xlabel('Reads Fraction')
    axes[1].set_ylabel('Reads')
    axes[1].set_yscale('log')
    plt.tight_layout()
    plt.savefig(f"{output}/reads_fraction_scatter.png")

    return df_major