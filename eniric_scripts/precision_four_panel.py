"""Plot a Figueira et al. 2016 like plot."""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import eniric
from eniric.utilities import rv_cumulative, rv_cumulative_full

matplotlib.use("Agg")


def load_dataframe(filename):
    """Load in precision file, clean up spaces in csv."""
    df = pd.read_csv(precision_file)
    # Temp, logg, [Fe/H], Alpha, Band, Resolution, vsini, Sampling, Quality, Cond. 1, Cond. 2, Cond. 3, correct flag
    df.columns = df.columns.str.strip()
    df.Band = df.Band.str.strip()
    df.Resolution = df.Resolution.str.strip()
    df.Quality = df.Quality.astype(float)
    df.Temp = df.Temp.astype(float)
    df.logg = df.logg.astype(float)
    df["[Fe/H]"] = df["[Fe/H]"].astype(float)
    df.Alhpa = df.Alpha.astype(float)
    df.vsini = df.vsini.astype(float)
    df.Sampling = df.Sampling.astype(float)
    df["Cond. 1"] = df["Cond. 1"].astype(float)
    df["Cond. 2"] = df["Cond. 2"].astype(float)
    df["Cond. 3"] = df["Cond. 3"].astype(float)
    df["correct flag"] = df["correct flag"].astype(bool)

    return df


def plot_precision(
    precision_file, teffs=None, logg=4.5, fe_h=0.0, vsini=1.0, sampling=3
):
    """Plot precision 4 panel with RV. precision."""
    if teffs is None:
        # Default teffs
        teffs = [3900, 3500, 2800, 2600]
    assert len(teffs) == 4
    df = load_dataframe(precision_file)
    filter_dict = {"logg": logg, "[Fe/H]": fe_h, "vsini": vsini, "Sampling": sampling}
    df = filter_df(
        df, filter_dict, drop_list=["Alpha", "[Fe/H]", "correct flag", "Quality"]
    )

    fig, axes = plt.subplots(2, 2)
    ax = axes.flatten()
    df_selected = df[df.Resolution.isin(["60k", "80k", "100k"])]
    df_selected = df_selected[df_selected.Temp.isin(teffs)]

    maximums = []
    minimums = []

    for ii, temp in enumerate(teffs):
        # This entry
        df_ii = df_selected[df_selected["Temp"] == temp]
        df_ii_60k = df_ii[df_ii["Resolution"].str.strip() == "60k"]
        df_ii_80k = df_ii[df_ii["Resolution"].str.strip() == "80k"]
        df_ii_100k = df_ii[df_ii["Resolution"].str.strip() == "100k"]

        df_ii_60k = df_ii_60k.set_index("Band")
        print("new index\n", df_ii_60k)
        df_ii_60k = df_ii_60k.reindex(["Z", "Y", "J", "H", "K"])
        print("reindexed\n", df_ii_60k)

        df_ii_80k = df_ii_80k.set_index("Band")
        df_ii_80k = df_ii_80k.reindex(["Z", "Y", "J", "H", "K"])
        df_ii_100k = df_ii_100k.set_index("Band")
        df_ii_100k = df_ii_100k.reindex(["Z", "Y", "J", "H", "K"])

        maximums.append(
            np.max(
                [
                    df_ii_60k[["Cond. 1", "Cond. 2", "Cond. 3"]].max(),
                    df_ii_80k[["Cond. 1", "Cond. 2", "Cond. 3"]].max(),
                    df_ii_100k[["Cond. 1", "Cond. 2", "Cond. 3"]].max(),
                ]
            )
        )
        minimums.append(
            np.min(
                [
                    df_ii_60k[["Cond. 1", "Cond. 2", "Cond. 3"]].min(),
                    df_ii_80k[["Cond. 1", "Cond. 2", "Cond. 3"]].min(),
                    df_ii_100k[["Cond. 1", "Cond. 2", "Cond. 3"]].min(),
                ]
            )
        )

        ax[ii].fill_between(
            df_ii_60k.index,
            df_ii_60k["Cond. 2"].values,
            df_ii_60k["Cond. 3"].values,
            color="b",
            alpha=0.2,
        )
        ax[ii].fill_between(
            df_ii_80k.index,
            df_ii_80k["Cond. 2"].values,
            df_ii_80k["Cond. 3"].values,
            color="g",
            alpha=0.2,
        )
        ax[ii].fill_between(
            df_ii_100k.index,
            df_ii_100k["Cond. 2"].values,
            df_ii_100k["Cond. 3"].values,
            color="r",
            alpha=0.2,
        )

        ax[ii].plot(
            df_ii_60k.index, df_ii_60k["Cond. 1"].values, color="b", linestyle="--"
        )  # lim
        ax[ii].plot(
            df_ii_80k.index, df_ii_80k["Cond. 1"].values, color="g", linestyle="--"
        )  # lim
        ax[ii].plot(
            df_ii_100k.index, df_ii_100k["Cond. 1"].values, color="r", linestyle="--"
        )  # lim

        ax[ii].scatter(
            df_ii_60k.index,
            df_ii_60k["Cond. 2"].values,
            marker="^",
            color="b",
            alpha=0.4,
        )
        ax[ii].scatter(
            df_ii_60k.index,
            df_ii_60k["Cond. 3"].values,
            marker="o",
            color="b",
            alpha=0.4,
        )

        ax[ii].scatter(
            df_ii_80k.index,
            df_ii_80k["Cond. 3"].values,
            marker="^",
            color="g",
            alpha=0.4,
        )
        ax[ii].scatter(
            df_ii_80k.index,
            df_ii_80k["Cond. 2"].values,
            marker="o",
            color="g",
            alpha=0.4,
        )

        ax[ii].scatter(
            df_ii_100k.index,
            df_ii_100k["Cond. 3"].values,
            marker="^",
            color="r",
            alpha=0.4,
        )
        ax[ii].scatter(
            df_ii_100k.index,
            df_ii_100k["Cond. 2"].values,
            marker="o",
            color="r",
            alpha=0.4,
        )

    # Set limits ticks and labels
    ymax = np.max(maximums)
    ymin = np.min(minimums)
    delta_y = ymax - ymin
    band_size = len(df_ii_60k.index)
    for jj in range(4):
        ax[jj].text(0, ymax, "{} K".format(teffs[jj]), size=14)
        ax[jj].set_ylim(ymin - 0.11 * delta_y, ymax + 0.11 * delta_y)
        ax[jj].set_xlim(-0.5, band_size - 0.5)
        ax[jj].tick_params(axis="both", which="major", labelsize=12)

        # ticks and labels
        if (jj == 2) or (ii == 3):
            ax[jj].set_xlabel("Bands", fontsize=12)
        if (jj == 1) or (jj == 3):
            ax[jj].set_yticklabels([])
        if (jj == 0) or (jj == 1):
            ax[jj].set_xticklabels([])

    fig.text(
        0.04,
        0.5,
        r"Precision [m/s]",
        ha="center",
        va="center",
        rotation="vertical",
        size=12,
    )

    fig.subplots_adjust(hspace=0, wspace=0, bottom=0.12, top=0.95, right=0.95)

    fig.savefig("plots/precision_logg{0}_feh_{1}.pdf".format(logg, fe_h))
    fig.savefig("plots/precision_logg{0}_feh_{1}.png".format(logg, fe_h), dpi=400)


def filter_df(df, filter_dict, drop_list=None):
    """Filter DataFrame by dictionary of key and values."""
    for key, val in filter_dict.items():
        df = df[df[key] == val]
    if drop_list is not None:
        df = df.drop(columns=drop_list)
    return df


def cumulative_df(df, full_cum=False):
    """Index by band, """
    bands = df.index
    assert all(bands == ["Z", "Y", "J", "H", "K"]), bands
    if full_cum:
        cum_bands = ["Z", "ZY", "ZYJ", "ZYJH", "ZYJHK", "YJHK", "JHK", "HK", "K"]

        cum_dict = {
            "Band": cum_bands,
            "Cond. 1": rv_cumulative_full(df["Cond. 1"]),
            "Cond. 2": rv_cumulative_full(df["Cond. 2"]),
            "Cond. 3": rv_cumulative_full(df["Cond. 3"]),
        }
    else:
        cum_bands = ["Z", "ZY", "ZYJ", "ZYJH", "ZYJHK"]
        cum_dict = {
            "Band": cum_bands,
            "Cond. 1": rv_cumulative(df["Cond. 1"], single=True),
            "Cond. 2": rv_cumulative(df["Cond. 2"], single=True),
            "Cond. 3": rv_cumulative(df["Cond. 3"], single=True),
        }
    cum_df = pd.DataFrame(cum_dict)
    cum_df = cum_df.set_index("Band")
    cum_df = cum_df.reindex(cum_bands)
    return cum_df


def cumulative_plot(
    precision_file,
    teffs=None,
    logg=4.5,
    fe_h=0.0,
    vsini=1.0,
    sampling=3,
    full_cum=False,
):
    """RV precision with cumulative bands.

    full_cum: bool
    Cumlative over entire range [ "Z","ZY", "ZYJ", "ZYJH", "ZYJHK","YJHK", "JHK","HK","K"]

    """
    if teffs is None:
        # Default values
        teffs = [3900, 3500, 2800, 2600]
    assert len(teffs) == 4

    df = load_dataframe(precision_file)
    filter_dict = {"logg": logg, "[Fe/H]": fe_h, "vsini": vsini, "Sampling": sampling}
    df = filter_df(
        df, filter_dict, drop_list=["Alpha", "[Fe/H]", "correct flag", "Quality"]
    )

    fig, axes = plt.subplots(2, 2)
    ax = axes.flatten()
    df_selected = df[df.Resolution.isin(["60k", "80k", "100k"])]
    df_selected = df_selected[df_selected.Temp.isin(teffs)]

    maximums = []
    minimums = []

    for ii, temp in enumerate(teffs):
        # This entry
        df_ii = df_selected[df_selected["Temp"] == temp]
        df_ii_60k = df_ii[df_ii["Resolution"].str.strip() == "60k"]
        df_ii_80k = df_ii[df_ii["Resolution"].str.strip() == "80k"]
        df_ii_100k = df_ii[df_ii["Resolution"].str.strip() == "100k"]

        df_ii_60k = df_ii_60k.set_index("Band")
        print("new index\n", df_ii_60k)
        df_ii_60k = df_ii_60k.reindex(["Z", "Y", "J", "H", "K"])
        print("reindexed\n", df_ii_60k)

        df_ii_80k = df_ii_80k.set_index("Band")
        df_ii_80k = df_ii_80k.reindex(["Z", "Y", "J", "H", "K"])
        df_ii_100k = df_ii_100k.set_index("Band")
        df_ii_100k = df_ii_100k.reindex(["Z", "Y", "J", "H", "K"])

        # Cumulative
        df_ii_60k = cumulative_df(df_ii_60k, full_cum=full_cum)
        df_ii_80k = cumulative_df(df_ii_80k, full_cum=full_cum)
        df_ii_100k = cumulative_df(df_ii_100k, full_cum=full_cum)

        maximums.append(
            np.max(
                [
                    df_ii_60k[["Cond. 1", "Cond. 2", "Cond. 3"]].max(),
                    df_ii_80k[["Cond. 1", "Cond. 2", "Cond. 3"]].max(),
                    df_ii_100k[["Cond. 1", "Cond. 2", "Cond. 3"]].max(),
                ]
            )
        )
        minimums.append(
            np.min(
                [
                    df_ii_60k[["Cond. 1", "Cond. 2", "Cond. 3"]].min(),
                    df_ii_80k[["Cond. 1", "Cond. 2", "Cond. 3"]].min(),
                    df_ii_100k[["Cond. 1", "Cond. 2", "Cond. 3"]].min(),
                ]
            )
        )

        ax[ii].fill_between(
            df_ii_60k.index,
            df_ii_60k["Cond. 2"].values,
            df_ii_60k["Cond. 3"].values,
            color="b",
            alpha=0.2,
        )
        ax[ii].fill_between(
            df_ii_80k.index,
            df_ii_80k["Cond. 2"].values,
            df_ii_80k["Cond. 3"].values,
            color="g",
            alpha=0.2,
        )
        ax[ii].fill_between(
            df_ii_100k.index,
            df_ii_100k["Cond. 2"].values,
            df_ii_100k["Cond. 3"].values,
            color="r",
            alpha=0.2,
        )

        ax[ii].plot(
            df_ii_60k.index, df_ii_60k["Cond. 1"].values, color="b", linestyle="--"
        )  # lim
        ax[ii].plot(
            df_ii_80k.index, df_ii_80k["Cond. 1"].values, color="g", linestyle="--"
        )  # lim
        ax[ii].plot(
            df_ii_100k.index, df_ii_100k["Cond. 1"].values, color="r", linestyle="--"
        )  # lim

        ax[ii].scatter(
            df_ii_60k.index,
            df_ii_60k["Cond. 2"].values,
            marker="^",
            color="b",
            alpha=0.4,
        )
        ax[ii].scatter(
            df_ii_60k.index,
            df_ii_60k["Cond. 3"].values,
            marker="o",
            color="b",
            alpha=0.4,
        )

        ax[ii].scatter(
            df_ii_80k.index,
            df_ii_80k["Cond. 3"].values,
            marker="^",
            color="g",
            alpha=0.4,
        )
        ax[ii].scatter(
            df_ii_80k.index,
            df_ii_80k["Cond. 2"].values,
            marker="o",
            color="g",
            alpha=0.4,
        )

        ax[ii].scatter(
            df_ii_100k.index,
            df_ii_100k["Cond. 3"].values,
            marker="^",
            color="r",
            alpha=0.4,
        )
        ax[ii].scatter(
            df_ii_100k.index,
            df_ii_100k["Cond. 2"].values,
            marker="o",
            color="r",
            alpha=0.4,
        )

        ax[ii].fill_between(
            df_ii_100k.index,
            df_ii_100k["Cond. 2"].values,
            df_ii_100k["Cond. 3"].values,
            color="r",
            alpha=0.2,
        )

    # Set limits ticks and labels
    ymax = np.min([np.max(maximums), 10])  # Set limit to 20 km/s
    ymin = np.min(minimums)
    delta_y = ymax - ymin
    band_size = len(df_ii_60k.index)

    for jj in range(4):
        ax[jj].text(0, ymax, "{} K".format(teffs[jj]), size=14)
        ax[jj].set_ylim(ymin - 0.1 * delta_y, ymax + 0.15 * delta_y)
        ax[jj].set_xlim(-0.5, band_size - 0.5)
        ax[jj].tick_params(axis="both", which="major", labelsize=12)

        # ticks and labels
        if (jj == 2) or (ii == 3):
            ax[jj].set_xlabel("Bands", fontsize=12)
            if full_cum:
                ax[jj].tick_params(axis="x", labelrotation=40)
            else:
                ax[jj].tick_params(axis="x", labelrotation=25)

        if (jj == 1) or (jj == 3):
            ax[jj].set_yticklabels([])
        if (jj == 0) or (jj == 1):
            ax[jj].set_xticklabels([])

    fig.text(
        0.04,
        0.5,
        r"Precision [m/s]",
        ha="center",
        va="center",
        rotation="vertical",
        size=12,
    )

    fig.subplots_adjust(hspace=0, wspace=0, bottom=0.17, top=0.95, right=0.95)

    fname = "plots/cummulative_precision_logg{0}_feh_{1}_{2}".format(
        logg, fe_h, full_cum
    )
    fig.savefig(fname + ".pdf")
    fig.savefig(fname + ".png", dpi=400)


if __name__ == "__main__":
    # ToDo: Add an argparse cli

    precision_file = os.path.join(eniric.paths["precision"], "precision_results.csv")
    plot_precision(precision_file, teffs=[3900, 3500, 2800, 2600])
    # plot_precision(precision_file, teffs=[3900, 3500, 2800, 2600], logg=4, fe_h=1)

    cumulative_plot(precision_file, teffs=[3900, 3500, 2800, 2600])
    cumulative_plot(precision_file, teffs=[3900, 3500, 2800, 2600], full_cum=True)
    sys.exit(0)
