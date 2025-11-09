
import os
import glob
import pandas as pd
from enrichment_analyzer import run_enrichment
from overlap_visualizer import plot_heatmap, plot_network
from report_generator import generate_html

DATA_DIR = 'data'
RES_DIR = 'results'
GMT_PATH = os.path.join(DATA_DIR, 'genesets.gmt')


def main():
    os.makedirs(RES_DIR, exist_ok=True)
    frames = []
    for f in sorted(glob.glob(os.path.join(DATA_DIR, '*_genes.txt'))):
        print(f'🧠 Analyzing {f} ...')
        df = run_enrichment(f, GMT_PATH)
        frames.append(df)
    if not frames:
        raise SystemExit('No gene lists found in data/*.txt')
    all_df = pd.concat(frames, ignore_index=True)
    # Save tables
    out_tbl = os.path.join(RES_DIR, 'enrichment_tables', 'all_enrichments.csv')
    os.makedirs(os.path.dirname(out_tbl), exist_ok=True)
    all_df.to_csv(out_tbl, index=False)
    # Plots
    plot_heatmap(all_df, os.path.join(RES_DIR, 'figures', 'overlap_heatmap.png'))
    plot_network(all_df, os.path.join(RES_DIR, 'figures', 'pathway_network.png'))
    # Report
    generate_html(all_df, os.path.join(RES_DIR, 'reports', 'pathway_report.html'))
    print('✅ Analysis complete. Results stored in ./results')


if __name__ == '__main__':
    main()
