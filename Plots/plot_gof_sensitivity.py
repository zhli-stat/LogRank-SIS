import matplotlib.pyplot as plt
import numpy as np

d_values = [15, 30, 50, 60]

data = {
    'Log-rank SIS': {
        'LogLik': [-357.263, -331.806, -295.865, -285.005],
        'ChiSq': [138.830, 189.743, 261.626, 283.345],
        'AIC': [744.526, 721.613, 689.730, 686.011],
        'BIC': [798.940, 826.814, 867.483, 896.413]
    },
    'ADD-SIS': {
        'LogLik': [-387.517, -369.296, -350.239, -327.011],
        'ChiSq': [78.321, 114.764, 152.878, 199.334],
        'AIC': [805.034, 796.592, 794.478, 768.022],
        'BIC': [859.449, 901.793, 964.976, 974.796]
    },
    'MV-SIS (M3)': {
        'LogLik': [-386.409, -362.292, -329.998, -322.471],
        'ChiSq': [80.538, 128.771, 193.359, 208.415],
        'AIC': [804.817, 786.585, 759.997, 762.941],
        'BIC': [862.859, 899.041, 941.378, 976.971]
    },
    'DC-SIS (M1)': {
        'LogLik': [-387.240, -360.189, -341.954, -328.795],
        'ChiSq': [78.876, 132.977, 169.447, 195.766],
        'AIC': [806.480, 782.379, 785.909, 777.590],
        'BIC': [864.522, 894.835, 970.918, 995.247]
    },
    'KF (M2)': {
        'LogLik': [-391.735, -363.075, -343.406, -329.948],
        'ChiSq': [69.885, 127.205, 166.544, 193.460],
        'AIC': [807.470, 782.150, 782.811, 773.895],
        'BIC': [851.002, 883.724, 956.937, 980.670]
    }
}

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