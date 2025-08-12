#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mutations", required=True, help="mutations_long.csv")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--min_count", type=int, default=20, help="filter mutations by total count")
    args = p.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.mutations, parse_dates=["collection_date"])
    df = df.dropna(subset=["collection_date"])

    # total sequences per day
    day_totals = df.groupby(df["collection_date"].dt.date)["key"].nunique().rename("n_sequences")
    day_totals = day_totals.reset_index().rename(columns={"collection_date":"date"})

    # mutation counts per day
    day_mut = df.groupby([df["collection_date"].dt.date, "label"])["key"].nunique().reset_index()
    day_mut = day_mut.rename(columns={"collection_date":"date", "key":"n_with_mut"})

    # merge to compute frequency
    merged = day_mut.merge(day_totals, on="date", how="left")
    merged["freq"] = merged["n_with_mut"] / merged["n_sequences"].replace({0: pd.NA})
    # filter rare
    top = merged.groupby("label")["n_with_mut"].sum().sort_values(ascending=False).reset_index()
    keep = top[top["n_with_mut"] >= args.min_count]["label"]
    merged = merged[merged["label"].isin(keep)]

    merged.sort_values(["label","date"], inplace=True)
    merged.to_csv(outdir / "mutation_freq_by_day.csv", index=False)
    print(f"Wrote {outdir/'mutation_freq_by_day.csv'}")

if __name__ == "__main__":
    main()
