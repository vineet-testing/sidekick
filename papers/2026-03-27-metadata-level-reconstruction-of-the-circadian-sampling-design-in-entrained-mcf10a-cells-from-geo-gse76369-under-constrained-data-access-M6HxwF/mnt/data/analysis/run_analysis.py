
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import requests
import textwrap

base = Path('/mnt/data')
(base/'logs').mkdir(exist_ok=True)
(base/'tables').mkdir(exist_ok=True)
(base/'analysis').mkdir(exist_ok=True)
(base/'figures').mkdir(exist_ok=True)

# Metadata reconstructed from the provided dataset profile
samples = [
    {"gsm":"GSM1982275","time_h":0,"gpr_file":"GSE76369_0042.gpr","channel":"Cy3","pair_time_h":16},
    {"gsm":"GSM1982276","time_h":4,"gpr_file":"GSE76369_0219.gpr","channel":"Cy3","pair_time_h":20},
    {"gsm":"GSM1982277","time_h":8,"gpr_file":"GSE76369_0041.gpr","channel":"Cy3","pair_time_h":24},
    {"gsm":"GSM1982278","time_h":12,"gpr_file":"GSE76369_0220.gpr","channel":"Cy3","pair_time_h":28},
    {"gsm":"GSM1982279","time_h":16,"gpr_file":"GSE76369_0042.gpr","channel":"Cy5","pair_time_h":0},
    {"gsm":"GSM1982280","time_h":20,"gpr_file":"GSE76369_0219.gpr","channel":"Cy5","pair_time_h":4},
    {"gsm":"GSM1982281","time_h":24,"gpr_file":"GSE76369_0041.gpr","channel":"Cy5","pair_time_h":8},
    {"gsm":"GSM1982282","time_h":28,"gpr_file":"GSE76369_0220.gpr","channel":"Cy5","pair_time_h":12},
]
df = pd.DataFrame(samples)
df["cell_line"] = "MCF10A"
df["phenotype_harmonized"] = "epithelial, non-tumorigenic / noncancerous"
df.to_csv(base/'tables'/'sample_metadata.tsv', sep='\t', index=False)

pair_rows = []
for gpr, sub in df.groupby('gpr_file'):
    cy3 = sub.loc[sub['channel'].eq('Cy3')].iloc[0]
    cy5 = sub.loc[sub['channel'].eq('Cy5')].iloc[0]
    pair_rows.append({
        'gpr_file': gpr,
        'cy3_gsm': cy3['gsm'],
        'cy3_time_h': int(cy3['time_h']),
        'cy5_gsm': cy5['gsm'],
        'cy5_time_h': int(cy5['time_h']),
        'time_offset_h': int(cy5['time_h'] - cy3['time_h']),
    })
pairs = pd.DataFrame(pair_rows).sort_values('cy3_time_h')
pairs.to_csv(base/'tables'/'raw_pairing_design.tsv', sep='\t', index=False)

# Attempt network retrieval required for the selected experiments.
targets = [
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76369",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21250",
] + [f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}" for gsm in df['gsm']]

results = []
for url in targets:
    rec = {"url": url}
    try:
        r = requests.get(url, timeout=20)
        rec["status_code"] = r.status_code
        rec["ok"] = bool(r.ok)
        rec["n_bytes"] = len(r.content)
    except Exception as e:
        rec["ok"] = False
        rec["error"] = repr(e)
    results.append(rec)
network_df = pd.DataFrame(results)
network_df.to_csv(base/'logs'/'network_attempts.tsv', sep='\t', index=False)

# Metadata-level timeline figure as the same-dataset pivot.
plt.figure(figsize=(9, 2.8))
plt.hlines(1, xmin=0, xmax=28, color='gray', linewidth=2)
color_map = {'Cy3': '#1f77b4', 'Cy5': '#d62728'}
plt.scatter(df['time_h'], [1]*len(df), c=[color_map[c] for c in df['channel']], s=100, zorder=3)
for _, row in df.iterrows():
    y = 1.08 if row['channel'] == 'Cy3' else 0.90
    va = 'bottom' if row['channel'] == 'Cy3' else 'top'
    plt.text(row['time_h'], y, f"{row['time_h']}h\n{row['gsm'][-4:]}", ha='center', va=va, fontsize=8)
for _, row in pairs.iterrows():
    plt.plot([row['cy3_time_h'], row['cy5_time_h']], [1, 1], color='black', alpha=0.3, linewidth=4, solid_capstyle='butt')
    midpoint = (row['cy3_time_h'] + row['cy5_time_h']) / 2
    label = row['gpr_file'].replace('GSE76369_', '').replace('.gpr', '')
    plt.text(midpoint, 1.18, label, ha='center', fontsize=8)
plt.yticks([])
plt.xticks(range(0, 29, 4))
plt.xlabel('Time after serum shock (h)')
plt.title('GSE76369 MCF10A sample schedule and two-color raw pairings')
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Cy3 sample', markerfacecolor=color_map['Cy3'], markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Cy5 sample', markerfacecolor=color_map['Cy5'], markersize=8),
]
plt.legend(handles=legend_elements, loc='upper left', frameon=False)
plt.ylim(0.72, 1.3)
plt.tight_layout()
plt.savefig(base/'figures'/'figure_1_design_timeline.png', dpi=160, bbox_inches='tight')
plt.close()

# Narrative handoff
summary = textwrap.dedent("""
# GSE76369 offline execution summary

## What was attempted
1. Re-download GEO resources needed for the selected experiments:
   - sample-level processed expression tables for GSM1982275-GSM1982282
   - GPL21250 platform annotation
   - series-level metadata pages
2. If successful, reconstruct the 44,544-probe by 8-sample matrix and run conservative temporal screening.

## Outcome
The compute sandbox has no internet access. Re-download attempts to NCBI GEO failed, and no cached GSE76369 data files were present locally in `/mnt/data`. Because the numeric processed tables were unavailable, the selected matrix-reconstruction and temporal-ranking experiments could not be completed faithfully.

## Pivot completed on the same dataset
A metadata-only handoff was produced from the provided dataset profile:
- `tables/sample_metadata.tsv`: timepoints, GSM accessions, raw-file names, and channel assignments
- `tables/raw_pairing_design.tsv`: four 16-hour two-color pairings
- `figures/figure_1_design_timeline.png`: visual overview of the sample schedule and raw pairings
- `logs/network_attempts.tsv`: reproducible record of failed re-download attempts

## Practical next step
Re-run `analysis/run_analysis.py` in a network-enabled environment or with the GEO sample tables and GPL21250 annotation mounted locally; then proceed with probe-matrix reconstruction, annotation QC, and descriptive 24-hour harmonic ranking.
""").strip()
(base/'analysis_summary.md').write_text(summary)
