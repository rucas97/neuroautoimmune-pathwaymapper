
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx


def plot_heatmap(all_df: pd.DataFrame, out_path: str):
    # Build matrix of -log10(FDR) per pathway x disease
    mat = all_df.pivot_table(index='Term', columns='disease', values='Adjusted P-value', aggfunc='min')
    mat = mat.applymap(lambda x: -np.log10(x) if isinstance(x, (int, float)) and x>0 else 0)
    plt.figure(figsize=(10, max(6, 0.3*len(mat))))
    sns.heatmap(mat.fillna(0), cmap='viridis')
    plt.title('Shared NeuroAutoImmune Pathway Enrichment (-log10 FDR)')
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=220)
    plt.close()


def plot_network(all_df: pd.DataFrame, out_path: str, fdr_thresh: float = 0.1):
    # Bipartite graph: disease <-> pathway edges for significant terms
    G = nx.Graph()
    sig = all_df[all_df['Adjusted P-value'] <= fdr_thresh]
    for _, r in sig.iterrows():
        d = r['disease']
        t = r['Term']
        w = max(0.1, -np.log10(max(r['Adjusted P-value'], 1e-300)))
        G.add_node(d, bipartite=0)
        G.add_node(t, bipartite=1)
        G.add_edge(d, t, weight=w)
    plt.figure(figsize=(12, 8))
    if G.number_of_nodes() == 0:
        plt.text(0.5,0.5,'No significant pathways at FDR<=%.2f' % fdr_thresh, ha='center')
    else:
        pos = nx.spring_layout(G, k=0.4, seed=42)
        nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', font_size=8)
    plt.title('NeuroAutoImmune Pathway Overlap Network (FDR<=%.2f)' % fdr_thresh)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=220)
    plt.close()
