
import os
import pandas as pd
from scipy.stats import hypergeom


def read_gmt(gmt_path: str) -> dict:
    genesets = {}
    with open(gmt_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('	')
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = set([g.strip() for g in parts[2:] if g.strip()])
            if genes:
                genesets[name] = genes
    return genesets


def bh_fdr(pvals):
    # Benjamini-Hochberg
    import numpy as np
    p = np.asarray(pvals)
    n = p.size
    order = np.argsort(p)
    ranks = np.empty(n, int)
    ranks[order] = np.arange(1, n+1)
    adj = p * n / ranks
    # ensure monotonicity
    adj_sorted = adj[order]
    for i in range(n-2, -1, -1):
        adj_sorted[i] = min(adj_sorted[i], adj_sorted[i+1])
    out = adj_sorted.copy()
    out[order] = adj_sorted  # map back
    # cap at 1
    out = out.clip(max=1.0)
    return out


def run_enrichment(gene_file: str, gmt_path: str) -> pd.DataFrame:
    # Load user genes
    genes = set(g.strip() for g in open(gene_file, 'r', encoding='utf-8') if g.strip())
    disease = os.path.basename(gene_file).replace('_genes.txt', '')

    # Load genesets
    gsets = read_gmt(gmt_path)
    # Universe = all genes appearing in any set union genes
    universe = set().union(*gsets.values()) | genes
    N = len(universe)
    n = len(genes)

    rows = []
    for term, term_genes in gsets.items():
        K = len(term_genes)
        k = len(genes & term_genes)
        if k == 0:
            pval = 1.0
        else:
            # P(X >= k) in hypergeom
            rv = hypergeom(N, K, n)
            pval = rv.sf(k-1)
        rows.append((term, K, k, pval, ','.join(sorted(genes & term_genes))))

    df = pd.DataFrame(rows, columns=['Term','SetSize','OverlapCount','P-value','OverlapGenes'])
    df['Adjusted P-value'] = bh_fdr(df['P-value'].values)
    df['disease'] = disease
    return df
