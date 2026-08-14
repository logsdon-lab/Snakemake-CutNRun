import os
import polars as pl
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from collections.abc import Sequence

KWARGS = {
    "separator": "\t",
    "columns": [0, 1],
    "has_header": False,
    "new_columns": ["chrom", "length"]
}

WD = os.path.dirname(__file__)

def fmt_df(old: str, new: str) -> pl.DataFrame:
    df_ctrl_new = pl.read_csv(new, **KWARGS)
    df_ctrl_old = pl.read_csv(old, **KWARGS)
    return (
        df_ctrl_old.join(df_ctrl_new, how="full", on="chrom")
        .drop("chrom_right")
        .rename({"length_right": "length_new"})
        .fill_null(0)
    )

def draw_hist(ax: Axes, srs: pl.Series, title: str) -> None:
    max_cnt = srs.value_counts()["count"].max()
    mode = srs.mode().to_list()
    if isinstance(max_cnt, (float, int)):
        ypos = max_cnt / 2
    else:
        ypos = 1000

    for v in mode:
        ax.axvline(v, color="black", label="mode", linestyle="dotted")
    median = srs.median()
    if median:
        assert isinstance(median, (float, int))
        
        ax.annotate(f"{median}", (median, ypos), color="red")
        ax.axvline(median, color="red", label="median", linestyle="dotted")

    print(f"{mode=}, {median=}")

    ax.set_title(title)
    ax.set_xlabel("Difference in length (bp)")
    ax.set_ylabel("Number of reads")
    ax.hist(srs, bins=100)

df_ctrl = fmt_df(
    f"{WD}/../../results/IgG/trimmed_reads/reads_trimmed_filtered.fq.fai",
    f"{WD}/../../results_new/NA20355/trimmed_reads/reads_control_trimmed.fq.gz.fai"
)
srs_ctrl_length_diff = df_ctrl.select(diff=pl.col("length_new")-pl.col("length"))["diff"]

df_trt = fmt_df(
    f"{WD}/../../results/CENP-A/trimmed_reads/reads_trimmed_filtered.fq.fai",
    f"{WD}/../../results_new/NA20355/trimmed_reads/reads_treatment_trimmed.fq.gz.fai",
)
srs_trt_length_diff = df_trt.select(diff=pl.col("length_new")-pl.col("length"))["diff"]

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, layout="constrained")
axes: Sequence[Axes]

draw_hist(axes[0], srs_ctrl_length_diff, "Control (IgG)")
axes[0].set_xlabel("")

draw_hist(axes[1], srs_trt_length_diff, "Treatment (CENP-A)")

h, l = axes[0].get_legend_handles_labels()
fig.legend(handles=h, labels=l)
fig.savefig(os.path.join(WD, "hist_read_length_cmp.png"), dpi=600, bbox_inches="tight")
