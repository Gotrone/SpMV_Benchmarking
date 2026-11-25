"""Plot the Task 2 timings (16384 x 16384 matrix, 10 % non-zeros) reported in Task2/README.md.

Run from the repository root:  python scripts/make_figures.py
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

levels = ["-O0", "-O2 no-vec", "-O3 vec", "-Ofast vec"]
kernels = ["Dense", "COO", "CSR", "CSC"]
# execution time in ms, copied from Task2/README.md
data = {
    "GCC": {"Dense": [764, 355, 366, 782], "COO": [97, 87, 72, 72], "CSR": [85, 74, 131, 74], "CSC": [105, 33, 30, 31]},
    "ICC": {"Dense": [749, 357, 157, 158], "COO": [97, 75, 73, 73], "CSR": [77, 72, 74, 183], "CSC": [92, 28, 31, 30]},
}
colors = ["#7f8c8d", "#2980b9", "#27ae60", "#c0392b"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
x = np.arange(len(levels))
w = 0.2
for ax, (compiler, vals) in zip(axes, data.items()):
    for i, k in enumerate(kernels):
        ax.bar(x + (i - 1.5) * w, vals[k], w, label=k, color=colors[i])
    ax.set_xticks(x, levels)
    ax.set_title(f"{compiler}")
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("Execution time (ms, lower is better)")
axes[1].legend(title="Kernel")
fig.suptitle("SpMV: dense vs COO / CSR / CSC, 16384 x 16384 matrix, 10 % non-zeros")
fig.tight_layout()
fig.savefig("images/spmv_timings.png", dpi=150)
