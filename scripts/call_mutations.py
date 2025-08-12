#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from datetime import datetime
from dateutil.parser import parse as dtparse
from pathlib import Path

def load_metadata(meta_path):
    meta = pd.read_csv(meta_path)
    # normalize date column
    date_col = None
    for cand in ["collection_date", "date", "Collection_Date", "sample_date"]:
        if cand in meta.columns:
            date_col = cand
            break
    if date_col is None:
        raise ValueError("Metadata must include a collection date column.")
    meta["collection_date"] = pd.to_datetime(meta[date_col].apply(dtparse), errors="coerce")
    # accession/id harmonization best-effort
    for cand in ["accession", "Accession", "id", "ID", "strain"]:
        if cand in meta.columns:
            meta["key"] = meta[cand].astype(str)
            break
    if "key" not in meta.columns:
        # fallback index as key
        meta["key"] = meta.index.astype(str)
    return meta[["key","collection_date"]]

def align_and_mutations(ref_seq, query_seq):
    # Simple global alignment (no gaps penalties tuned)
    alns = pairwise2.align.globalms(ref_seq, query_seq, 2, -1, -2, -0.5, one_alignment_only=True)
    aln_ref, aln_qry, score, start, end = alns[0]
    mutations = []
    ref_pos = 0
    for i, (r, q) in enumerate(zip(aln_ref, aln_qry)):
        if r != "-":
            ref_pos += 1
        if r == "-" or q == "-":
            # skip indels for frequency tracking (optional: include as 'del/ins')
            continue
        if r != q:
            mutations.append((ref_pos, r, q))  # 1-based ref position
    return mutations

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True, help="Path to sequences FASTA")
    p.add_argument("--metadata", required=True, help="CSV with accession + collection_date")
    p.add_argument("--out", required=True, help="Output folder")
    args = p.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(args.fasta, "fasta"))
    if len(records) == 0:
        raise ValueError("No sequences found in FASTA.")
    ref = str(records[0].seq)  # take first as reference (customize as needed)

    meta = load_metadata(args.metadata)

    rows = []
    for rec in records:
        key = rec.id
        seq = str(rec.seq)
        muts = align_and_mutations(ref, seq)
        # date lookup
        hit = meta.loc[meta["key"] == key]
        date = pd.NaT if hit.empty else hit["collection_date"].iloc[0]
        for (pos, ref_base, alt_base) in muts:
            rows.append({
                "key": key,
                "collection_date": date,
                "pos": pos,
                "ref": ref_base,
                "alt": alt_base,
                "label": f"{ref_base}{pos}{alt_base}"
            })

    mut_long = pd.DataFrame(rows)
    mut_long.to_csv(outdir / "mutations_long.csv", index=False)
    print(f"Wrote {outdir/'mutations_long.csv'} with {len(mut_long)} rows")

if __name__ == "__main__":
    main()
