import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path
import pandas as pd

IN_CSV = (Path(__file__).resolve().parents[1]
          / "CaseStudy" / "Output" / "ADNI" / "adni_gof_by_d.csv")
if not IN_CSV.exists():
    raise SystemExit(f"Input not found: {IN_CSV}. Knit the ADNI Rmd first.")
gof = pd.read_csv(IN_CSV)

SHOWN = {
    'Log-rank SIS': 'Log-rank SIS',
    'ADD-SIS':      'ADD-SIS',
    'MV-SIS (M3)':  'MV-SIS-M3',
    'DC-SIS (M1)':  'DC-SIS-M1',
    'KF (M2)':      'KS-Filter-M2',
}

d_values = sorted(gof['d'].unique().tolist())
data = {}
for label, method in SHOWN.items():
    sub = gof[gof['Method'] == method].set_index('d').reindex(d_values)
    data[label] = {m: sub[m].tolist() for m in ['LogLik', 'ChiSq', 'AIC', 'BIC']}

print(gof[gof['d'] == 50][['Method', 'LogLik']].to_string(index=False))

colors = {
    'Log-rank SIS': '#c44e52',
    'ADD-SIS': '#4c72b0',
    'MV-SIS (M3)': '#55a868',
    'DC-SIS (M1)': '#dd8452',
    'KF (M2)': '#8172b3'
}

linestyles = {
    'Log-rank SIS': '-',
    'ADD-SIS': '--',
    'MV-SIS (M3)': '-.',
    'DC-SIS (M1)': ':',
    'KF (M2)': (0, (5, 2))
}

markers = {
    'Log-rank SIS': 'o',
    'ADD-SIS': 'v',
    'MV-SIS (M3)': '*',
    'DC-SIS (M1)': 's',
    'KF (M2)': '^'
}

marker_sizes = {
    'Log-rank SIS': 8,
    'ADD-SIS': 7,
    'MV-SIS (M3)': 11,
    'DC-SIS (M1)': 7,
    'KF (M2)': 7
}

metrics = ['LogLik', 'ChiSq', 'AIC', 'BIC']
y_labels = [
    r'Log-likelihood ($\uparrow$)',
    r'$\chi^2$ Statistic ($\uparrow$)',
    r'AIC ($\downarrow$)',
    r'BIC ($\downarrow$)'
]
panel_labels = ['(a)', '(b)', '(c)', '(d)']

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.weight': 'normal',
    'axes.labelweight': 'normal',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
})

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
axes = axes.flatten()

lines, legend_labels = [], []

for i, (metric, ylabel, plabel) in enumerate(
        zip(metrics, y_labels, panel_labels)):
    ax = axes[i]

    for method in data:
        lw = 2.5 if method == 'Log-rank SIS' else 1.8
        zorder = 10 if method == 'Log-rank SIS' else 5

        line, = ax.plot(
            d_values,
            data[method][metric],
            color=colors[method],
            linestyle=linestyles[method],
            marker=markers[method],
            linewidth=lw,
            markersize=marker_sizes[method],
            zorder=zorder
        )

        if i == 0:
            lines.append(line)
            legend_labels.append(method)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.yaxis.grid(True, linestyle='--', linewidth=0.5,
                  color='gray', alpha=0.35)
    ax.set_axisbelow(True)

    ax.set_xticks(d_values)
    ax.set_xlabel('Submodel Size $d$', fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(axis='both', which='major',
                   labelsize=12, direction='out', length=4)

    ax.text(0.03, 0.97, plabel,
            transform=ax.transAxes,
            fontsize=13, va='top', ha='left')

fig.legend(
    lines, legend_labels,
    loc='lower center',
    bbox_to_anchor=(0.5, 0.0),
    ncol=3,
    prop={'size': 12},
    frameon=False,
    columnspacing=2.0,
    handletextpad=0.5,
    handlelength=3.0
)

plt.tight_layout(rect=[0, 0.12, 1, 1], w_pad=1.0, h_pad=1.0)

plt.savefig('GoF_Sensitivity_Analysis.pdf',
            format='pdf', bbox_inches='tight')
plt.savefig('GoF_Sensitivity_Analysis.png',
            dpi=600, bbox_inches='tight')
plt.show()