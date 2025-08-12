# COVID-19 Mutation Trend Tracker

Track SARS-CoV-2 single-nucleotide mutation frequencies over time from public genome datasets.

## Dataset
Use a Kaggle dataset with sequences + metadata (date). Example: 
- https://www.kaggle.com/datasets/pranavraikokte/covid19-genome-sequences

Place the FASTA/CSV and metadata in `data/raw/`:
- `sequences.fasta` (or `sequences.csv` with a `sequence` column)
- `metadata.csv` with at least `accession` and `collection_date` (YYYY-MM-DD)

## What this repo includes
- **Notebook**: `notebooks/01_parse_and_trend.ipynb` for EDA, mutation calling vs reference, and time-series plots
- **Scripts**: 
  - `scripts/call_mutations.py` — Parse sequences, align to a reference, compute SNPs, write per-sample mutation tables
  - `scripts/build_trends.py` — Aggregate per-day mutation frequencies and save CSVs + figures
- **App**: `app/streamlit_app.py` — Minimal Streamlit app to explore mutation frequencies interactively
- **Results**: Figures and processed CSVs saved under `results/` and `data/processed/`

## Quickstart
```bash
# 1) Install deps
pip install -r requirements.txt

# 2) Put your files in data/raw/
#    sequences.fasta, metadata.csv

# 3) Run mutation calling
python scripts/call_mutations.py --fasta data/raw/sequences.fasta --metadata data/raw/metadata.csv --out data/processed

# 4) Build frequency trends
python scripts/build_trends.py --mutations data/processed/mutations_long.csv --out results

# 5) (Optional) Launch the app
streamlit run app/streamlit_app.py
```

## Notes
- Alignment uses a simple global alignment (Biopython pairwise2 with match/mismatch scoring). For large datasets, consider faster tools (minimap2 + maf parsing).
- The notebook uses **matplotlib** for all plots.
- Dates are parsed as daily; you can resample to weekly/monthly in the notebook.
