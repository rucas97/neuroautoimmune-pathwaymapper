
# 🧠 NeuroAutoImmune‑PathwayMapper

Author: Reza Anvaripour (rucas97)
License: MIT

A reproducible Python + Docker pipeline detecting shared pathways among autoimmune and neurodegenerative diseases using a local pathway library and over‑representation analysis (hypergeometric, BH‑FDR).

## 🚀 Quick Start
```bash
docker build -t neuro-pathmap .
docker run --rm -v ${PWD}:/app neuro-pathmap
```
Outputs go to `results/`:
- `enrichment_tables/all_enrichments.csv`
- `figures/overlap_heatmap.png`
- `figures/pathway_network.png`
- `reports/pathway_report.html`

## 🧬 Notes
- No internet required. Uses bundled `data/genesets.gmt`.
- ORA via scipy hypergeometric; FDR via Benjamini–Hochberg.

## 📂 Structure
```
data/       # input gene lists + genesets.gmt
src/        # analysis scripts
results/    # auto-generated outputs
```

## 📚 Citation
Anvaripour, R. (2025). NeuroAutoImmune‑PathwayMapper. GitHub: https://github.com/rucas97
