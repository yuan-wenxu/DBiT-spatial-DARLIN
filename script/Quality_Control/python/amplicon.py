from preprocessing import open_fastq_file, iter_fastq_raw
from utils import group_fraction, k_disperse
from tqdm import tqdm
from umi_tools import UMIClusterer
import pandas as pd
import matplotlib.pyplot as plt
import math
import argparse

class CorrectionConfig:
    def __init__(self, umi_hd_threshold, lb_error_rate, major_fraction_threshold_molecule, reads_cutoff, slope_cutoff):
        self.umi_hd_threshold = umi_hd_threshold                     # edit-distance threshold for UMI clustering within each SR
        self.lb_error_rate = lb_error_rate                     # per-base error rate used to derive LB edit-distance threshold
        self.major_fraction_threshold_molecule = major_fraction_threshold_molecule  # per (CR, UR) group: keep LR with reads fraction >= this value
        self.reads_cutoff = reads_cutoff                        # only keep (SR, UR, LR) with supported reads >= this value
        self.slope_cutoff = slope_cutoff                        # only keep SR (spots) with k = reads/UMIs >= this value

def str_to_bool(value):
    """Convert str to bool"""
    if isinstance(value, bool):
        return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError(f'Boolean value expected, got: {value}')

def main(darlin_reads, bc_umi_reads, darlin, output_path, config):

    # Read DARLIN sequences
    if darlin:
        LB_list = []
        with open_fastq_file(darlin_reads) as handle:
            for _, seq, _ in tqdm(iter_fastq_raw(handle), desc="loading DARLIN fastq", unit_scale=True, unit=" reads"):
                LB_list.append(seq)

    # Read spatial barcode and UMI sequences
    SR_list, UB_list = [], []
    with open_fastq_file(bc_umi_reads) as handle:
        for _, seq, _ in tqdm(iter_fastq_raw(handle), desc="loading BC_UMI fastq", unit_scale=True, unit=" reads"):
            SR = seq[:16]
            ub = seq[16:]
            SR_list.append(SR)
            UB_list.append(ub)

    # Build dataframe
    if darlin:
        df = pd.DataFrame({"SR": SR_list, "UB": UB_list, "LB": LB_list})
        df["LB_len"] = df["LB"].str.len()
    else:
        df = pd.DataFrame({"SR": SR_list, "UB": UB_list})

    if darlin:
        # Plot DARLIN sequence length distribution
        plt.figure(figsize=(4, 2))
        plt.hist(df["LB_len"], bins=range(1, 300, 1), edgecolor="black")
        plt.xlabel("DARLIN sequence length")
        plt.ylabel("Number of reads")
        plt.title("DARLIN sequence length distribution")
        plt.tight_layout()
        plt.savefig(f"{output_path}/DARLIN_seq_len_distribution.png")

    columns = df.columns.tolist()

    df_uniq = (
        df
        .groupby(columns)
        .size()
        .rename("reads")
        .reset_index()
        )
    
    total_reads = df_uniq["reads"].sum()
    print("Number of reads:", df_uniq["reads"].sum())
    print("Unique records:", len(df_uniq))
    print("Number of SR:", df_uniq["SR"].nunique())
    print("Number of UB:", df_uniq["UB"].nunique())
    if darlin:
        print("Number of LB:", df_uniq["LB"].nunique())
    print("\n")

    # UMI correction with umi_tools
    print("Using umi_tools for UMI correction...")
    clusterer = UMIClusterer(cluster_method="directional")
    umi_mapping = {}
    for sr_id, sub_df in tqdm(df_uniq.groupby("SR"), desc="Correcting UMIs with umi_tools"):
        if sr_id == "not assigned":
            continue
        # Build counts dictionary for this cell
        umi_counts_str = dict(zip(sub_df["UB"], sub_df["reads"]))
        if not umi_counts_str:
            continue
        umi_counts_bytes = {umi.encode(): count for umi, count in umi_counts_str.items()}
        umi_groups = clusterer(umi_counts_bytes, threshold = config.umi_hd_threshold)
        for group in umi_groups:
            representative_bytes = max(group, key=lambda u: umi_counts_bytes.get(u, 0))
            representative_str = representative_bytes.decode()
            for umi_bytes in group:
                umi_str = umi_bytes.decode()
                umi_mapping[(sr_id, umi_str)] = representative_str
    # Apply mapping (fallback to original UB if not clustered)
    df_uniq["UR"] = df_uniq.apply(
        lambda row: umi_mapping.get((row["SR"], row["UB"]), row["UB"]),
        axis=1,
        )
    print("\n")

    # LB correction with umi_tools
    if darlin:
        print("Using umi_tools for LB correction...")
        lb_mapping = {}
        for (sr_id, umi, lb_len), sub_df in tqdm(
            df_uniq.groupby(["SR", "UR", "LB_len"]),
            desc="Correcting LB with umi_tools",
            ):
            lb_counts_str = dict(zip(sub_df["LB"], sub_df["reads"]))
            if not lb_counts_str:
                continue
            lb_counts_bytes = {lb.encode(): count for lb, count in lb_counts_str.items()}
            len_aware_threshold = max(1, int(math.ceil(lb_len * config.lb_error_rate)))
            lb_groups = clusterer(lb_counts_bytes, threshold=len_aware_threshold)
            for group in lb_groups:
                representative_bytes = max(group, key=lambda u: lb_counts_bytes.get(u, 0))
                representative_str = representative_bytes.decode()
                for lb_bytes in group:
                    lb_str = lb_bytes.decode()
                    lb_mapping[(sr_id, umi, lb_len, lb_str)] = representative_str
        # Apply mapping (fallback to original LB if not clustered)
        df_uniq["LR"] = df_uniq.apply(
            lambda row: lb_mapping.get((row["SR"], row["UR"], row["LB_len"], row["LB"]), row["LB"]),
            axis=1,
            )
        print("\n")

    # Collapse to unique molecules
    if darlin:
        df_uniq = (
            df_uniq
            .groupby(["SR", "UR", "LR", "LB_len"])
            .agg({"reads": "sum"})
            .reset_index()
        )
    else:
        df_uniq = (  
            df_uniq
            .groupby(["SR", "UR"])
            .agg({"reads": "sum"})
            .reset_index()
        )

    # Major LR filtration per (SR, UR)
    if darlin:
        df_uniq = group_fraction(df_uniq, output_path, config)

    # k
    df_summary = k_disperse(df_uniq, output_path)

    df_final = df_uniq.merge(df_summary[["SR", "k"]], on="SR", how="left")
    df_final = df_final[(df_final["k"] >= config.slope_cutoff) & (df_final["reads"] >= config.reads_cutoff)]

    number_of_reads = df_final["reads"].sum()
    number_of_molecules = len(df_final)
    number_of_spots = df_final["SR"].nunique()
    if darlin:
        number_of_LBs = df_final["LR"].nunique()
    reads_of_COCA = total_reads - number_of_reads

    print("After low quality SR removal:")
    print(f'Number of reads: {number_of_reads:,}')
    print(f'Number of molecules: {number_of_molecules:,}')
    print(f'Number of spots: {number_of_spots:,}')
    if darlin:
        print(f'Number of lineage barcodes: {number_of_LBs:,}')
    print(f"Reads of capture-oligo carryover artifacts: {reads_of_COCA:,}")
    print("\n")

    if darlin:
        df_final["n_LR"] = df_final.groupby("SR")["LR"].transform("nunique")
        df_plot = df_final[['SR', 'n_LR']].drop_duplicates()
        plt.figure(figsize=(3, 2))
        plt.hist(df_plot["n_LR"], bins=range(1, 8), edgecolor="white")
        plt.xlabel("Number of LRs per SR")
        plt.ylabel("Number of SRs")
        plt.yscale("log")
        plt.tight_layout()
        plt.savefig(f"{output_path}/LR_per_SR.png")

    df_final.to_csv(f"{output_path}/final.csv", index=False)

    if darlin:
        df_features = pd.DataFrame({
            "id": df_final["LR"].drop_duplicates().values,
            "name": df_final["LR"].drop_duplicates().values,
            })
        df_features["type"] = "Lineage Barcode"
        df_features.to_csv(f"{output_path}/features.csv", index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="DARLIN correction pipeline")
    parser.add_argument("-dr", "--darlin_reads", type=str, help="path to DARLIN fastq file")
    parser.add_argument("-bu", "--bc_umi_reads", type=str, help="path to BC_UMI fastq file")
    parser.add_argument("-d", "--darlin", type=str_to_bool, help="whether DARLIN sequences are provided")
    parser.add_argument("-o", "--output_path", type=str, help="path to output directory")
    
    parser.add_argument("--umi_hd_threshold", type=int, default=1, help="edit-distance threshold for UMI clustering within each SR")
    parser.add_argument("--lb_error_rate", type=float, default=0.02, help="per-base error rate used to derive LB edit-distance threshold")
    parser.add_argument("--major_fraction_threshold_molecule", type=float, default=0.8, help="per (CR, UR) group: keep LR with reads fraction >= this value")
    parser.add_argument("--reads_cutoff", type=int, default=10, help="only keep (SR, UR, LR) with supported reads >= this value")
    parser.add_argument("--slope_cutoff", type=int, default=10, help="only keep SR (spots) with k = reads/UMIs >= this value")
    args = parser.parse_args()

    darlin_reads = args.darlin_reads
    bc_umi_reads = args.bc_umi_reads
    darlin = args.darlin
    output_path = args.output_path
    config = CorrectionConfig(
        args.umi_hd_threshold,
        args.lb_error_rate,
        args.major_fraction_threshold_molecule,
        args.reads_cutoff,
        args.slope_cutoff,
    )

    main(darlin_reads, bc_umi_reads, darlin, output_path, config)