from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from darlin_core import analyze_sequences


LABEL_INFO = {
    "CA": {
        "bank_file": "allele_bank_Gr_CA.csv.gz",
        "config": "Col1a1",
    },
    "RA": {
        "bank_file": "allele_bank_Gr_RA.csv.gz",
        "config": "Rosa",
    },
    "TA": {
        "bank_file": "allele_bank_Gr_TA.csv.gz",
        "config": "Tigre",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter DARLIN cellfiltered.csv files by allele bank only."
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing CA/RA/TA subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for outputs (default: input directory).",
    )
    parser.add_argument(
        "--bank-dir",
        help="Directory containing allele_bank_Gr_*.csv.gz files.",
    )
    parser.add_argument(
        "--min-sequence-length",
        type=int,
        default=20,
        help="Minimum sequence length passed to analyze_sequences.",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        choices=sorted(LABEL_INFO.keys()),
        default=sorted(LABEL_INFO.keys()),
        help="Labels to process.",
    )
    return parser.parse_args()


def load_bank_sequences(bank_path: Path) -> set[str]:
    if not bank_path.exists():
        raise SystemExit(f"Bank file not found: {bank_path}")
    
    if bank_path.suffix == ".gz":
        bank_df = pd.read_csv(bank_path, sep=",", dtype=str, keep_default_na=False, compression="gzip")
    else:
        bank_df = pd.read_csv(bank_path, sep=",", dtype=str, keep_default_na=False)
    if bank_df.empty:
        return set()

    if "allele" not in bank_df.columns:
        raise SystemExit(
            f"Bank file {bank_path} is missing required column: allele"
        )

    return {
        str(allele).strip()
        for allele in bank_df["allele"].tolist()
        if str(allele).strip()
    }


def _ensure_analysis_columns(analyze_result: pd.DataFrame, cache_file: Path) -> None:
    required_columns = {"query", "mutations"}
    missing_columns = required_columns - set(analyze_result.columns)
    if missing_columns:
        raise SystemExit(
            f"Analyze result {cache_file} is missing required columns: {missing_columns}"
        )


def load_or_analyze_sequences(
    sequences: list[str],
    cache_file: Path,
    config: str,
    min_sequence_length: int,
    log_prefix: str = "",
) -> pd.DataFrame:
    if not sequences:
        return pd.DataFrame(columns=["query", "mutations"])

    prefix = f"{log_prefix}: " if log_prefix else ""
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    cached_rows = pd.DataFrame()
    missing_sequences = sequences

    if cache_file.exists():
        try:
            cache_df = pd.read_csv(cache_file, dtype=str, keep_default_na=False)
            _ensure_analysis_columns(cache_df, cache_file)
            cached_rows = cache_df[cache_df["query"].isin(sequences)].copy()
            cached_queries = set(cached_rows["query"].tolist())
            missing_sequences = [seq for seq in sequences if seq not in cached_queries]
            print(
                f"{prefix}Using cached analyze_sequences results "
                f"({len(cached_rows)} sequences)"
            )
        except Exception as exc:
            print(f"{prefix}Cache read failed ({exc}), reanalyzing...")
            cached_rows = pd.DataFrame()
            missing_sequences = sequences

    analyzed_rows = pd.DataFrame(columns=["query", "mutations"])
    if missing_sequences:
        print(
            f"{prefix}Running analyze_sequences on {len(missing_sequences)} unique sequences..."
        )
        try:
            result_obj = analyze_sequences(
                missing_sequences,
                config=config,
                min_sequence_length=min_sequence_length,
                verbose=False,
            )
            analyzed_rows = result_obj.to_df()
            _ensure_analysis_columns(analyzed_rows, cache_file)
        except Exception as exc:
            print(f"{prefix}analyze_sequences failed ({exc}), returning empty results.")
            analyzed_rows = pd.DataFrame(columns=["query", "mutations"])

    combined = pd.concat([cached_rows, analyzed_rows], ignore_index=True)
    if combined.empty:
        return combined

    combined = combined.drop_duplicates(subset=["query"], keep="last").copy()

    if not cached_rows.empty or not analyzed_rows.empty:
        updated_cache = combined
        if cache_file.exists():
            try:
                full_cache = pd.read_csv(cache_file, dtype=str, keep_default_na=False)
                _ensure_analysis_columns(full_cache, cache_file)
                retained_cache = full_cache[~full_cache["query"].isin(sequences)]
                updated_cache = pd.concat(
                    [retained_cache, combined],
                    ignore_index=True,
                )
                updated_cache = updated_cache.drop_duplicates(
                    subset=["query"],
                    keep="last",
                )
            except Exception:
                updated_cache = combined

        updated_cache.to_csv(cache_file, index=False)
        print(f"{prefix}Cached results to {cache_file}")

    return combined


def filter_dataframe_by_allele_bank(
    frame: pd.DataFrame,
    bank_file: Path,
    cache_file: Path,
    config: str,
    min_sequence_length: int,
    lr_col: str = "LR",
    log_prefix: str = "",
) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()

    sequences = sorted(
        {
            str(sequence).strip()
            for sequence in frame[lr_col].tolist()
            if str(sequence).strip()
        }
    )
    if not sequences:
        return frame.iloc[0:0].copy()

    bank = load_bank_sequences(bank_file)
    analyze_result = load_or_analyze_sequences(
        sequences,
        cache_file,
        config=config,
        min_sequence_length=min_sequence_length,
        log_prefix=log_prefix,
    )
    if analyze_result.empty:
        return frame.iloc[0:0].copy()

    analyze_result = analyze_result.copy()
    analyze_result["out_of_bank"] = ~analyze_result["mutations"].isin(bank)
    analyze_result = analyze_result[analyze_result["out_of_bank"]].copy()
    analyze_result[lr_col] = analyze_result["query"]

    return frame.merge(
        analyze_result[[lr_col, "mutations"]],
        on=lr_col,
        how="inner",
    ).copy()


def process_label(
    input_dir: Path,
    output_dir: Path,
    label: str,
    bank_file: Path,
    config: str,
    min_sequence_length: int,
) -> None:
    input_path = input_dir / label / "cellfiltered.n_LR_le_count.csv"
    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")

    frame = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    if frame.empty:
        print(f"{label}: No data")
        return

    if "LR" not in frame.columns:
        raise SystemExit(f"Input file {input_path} is missing required column: LR")

    cache_dir = output_dir / label
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / ".analyzed_cache.csv"

    filtered = filter_dataframe_by_allele_bank(
        frame,
        bank_file=bank_file,
        cache_file=cache_file,
        config=config,
        min_sequence_length=min_sequence_length,
        lr_col="LR",
        log_prefix=label,
    )

    print(f"{label}: {len(frame)} rows -> {len(filtered)} rows after bank filtering")
    out_path = output_dir / label / "cellfiltered.n_LR_le_count.bank_filtered.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(out_path, index=False)
    print(f"  Saved: {out_path}")


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    bank_dir = Path(args.bank_dir)

    for label in args.labels:
        info = LABEL_INFO[label]
        print(f"\n=== Processing {label} ===")
        process_label(
            input_dir,
            output_dir,
            label,
            bank_file=bank_dir / info["bank_file"],
            config=info["config"],
            min_sequence_length=args.min_sequence_length,
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
