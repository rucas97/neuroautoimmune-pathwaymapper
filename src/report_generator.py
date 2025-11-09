
import os
import pandas as pd
from jinja2 import Template


def generate_html(all_df: pd.DataFrame, output_file: str, fdr_thresh: float = 0.1):
    # Prepare summary: per disease significant terms
    grouped = {
        d: df.sort_values('Adjusted P-value').reset_index(drop=True)
        for d, df in all_df.groupby('disease')
    }
    template_str = """
    <html>
    <head>
      <meta charset="utf-8" />
      <title>NeuroAutoImmune Pathway Report</title>
      <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; }
        th, td { border: 1px solid #ccc; padding: 4px 8px; }
      </style>
    </head>
    <body>
      <h2>Shared Pathway Analysis Report</h2>
      <p>Significance threshold: FDR ≤ {{ fdr_thresh }}</p>
      {% for dis, df in grouped.items() %}
        <h3>{{ dis }}</h3>
        <table>
          <tr><th>Term</th><th>SetSize</th><th>Overlap</th><th>P-value</th><th>FDR</th><th>OverlapGenes</th></tr>
          {% for r in df.itertuples() %}
            <tr>
              <td>{{ r.Term }}</td>
              <td>{{ r.SetSize }}</td>
              <td>{{ r.OverlapCount }}</td>
              <td>{{ '%.3e' % r._4 }}</td>
              <td>{{ '%.3e' % r._6 }}</td>
              <td>{{ r.OverlapGenes }}</td>
            </tr>
          {% endfor %}
        </table>
      {% endfor %}
    </body>
    </html>
    """
    html = Template(template_str).render(grouped=grouped, fdr_thresh=fdr_thresh)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)
