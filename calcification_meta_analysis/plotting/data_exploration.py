# general
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

from calcification_meta_analysis.plotting import (
    climatology as climatology_plotting,
)
from calcification_meta_analysis.plotting import (
    plot_config,
    plot_utils,
)
from calcification_meta_analysis.processing import climatology as climatology_processing


def plot_study_timeseries(
    df: pd.DataFrame, ax=None, colorby="core_grouping", dpi: int = 300
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Plot the temporal distribution of studies and observation counts.

    Args:
        df (pd.DataFrame): The dataframe containing the data.
        ax (matplotlib.axes.Axes): The axes to plot the study timeseries on.
        colorby (str): The column to color the study timeseries by.
        dpi (int): The dpi of the figure.

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: The figure and axes.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(9, 3), dpi=dpi)

    df.dropna(subset=["year"], inplace=True)

    # count occurrences of year of each doi
    year_counts = df.groupby("year")["doi"].nunique()

    small_fontsize = 0.8 * plt.rcParams["font.size"]
    ax.bar(
        year_counts.index,
        year_counts.values,
        color="royalblue",
        width=150,
        alpha=0.5,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_ylabel("Number of studies", fontsize=plt.rcParams["font.size"])
    ax.tick_params(axis="both", which="major", labelsize=small_fontsize)

    # set y-ticks to appear every 5 units
    max_count = year_counts.max()
    y_ticks = np.arange(0, max_count + 5, 5)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticks, fontsize=small_fontsize)
    ax.grid(axis="y", linestyle="--", alpha=0.6)
    # group df by selected group
    grouped_df = df.groupby([colorby, "year"])
    # sum total n for each unique group and by year
    group_counts = grouped_df["n"].sum()
    # plot number of observations each year, by selected group
    unique_group = df[colorby].unique()
    n_ax = ax.twinx()
    group_palette = sns.color_palette("colorblind", len(unique_group))
    group_color_map = dict(zip(unique_group, group_palette))

    for group in unique_group:
        group_data = group_counts[group_counts.index.get_level_values(colorby) == group]
        n_ax.scatter(
            group_data.index.get_level_values("year"),
            group_data.values,
            color=group_color_map[group],
            alpha=1,
            s=20,
            label=group,
            edgecolor="white",
            linewidth=0.5,
        )

    n_ax.set_ylabel(
        "Number of observations",
        rotation=-90,
        labelpad=12,
        fontsize=plt.rcParams["font.size"],
    )
    n_ax.tick_params(axis="y", labelsize=small_fontsize)

    # align y-axis ticks with the bar chart
    n_yticks = np.arange(0, 4001, 1000)  # manual assignment of ticks
    n_ax.set_yticks(n_yticks)
    n_ax.set_yticklabels(n_yticks, fontsize=small_fontsize)

    legend = n_ax.legend(
        title=colorby.capitalize().replace("_", " "),
        loc="upper left",
        fontsize=small_fontsize,
        framealpha=0.7,
        title_fontsize=small_fontsize,
    )
    legend.get_title().set_ha("left")
    plt.xlabel("Year", fontsize=plt.rcParams["font.size"])
    plt.tight_layout()
    return fig, ax


def plot_effect_size_distributions(
    df,
    effect_sizes=[
        "cohens_d",
        "hedges_g",
        "relative_calcification",
        "absolute_calcification",
        "st_relative_calcification",
        "st_absolute_calcification",
    ],
    outlier_limits=None,
    title="Effect Size Distributions",
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Plots histograms and boxplots for the given effect sizes in the dataframe.

    Args:
        df (pd.DataFrame): The dataframe containing the effect size data.
        effect_sizes (list): List of effect size column names to plot.
        outlier_limits (tuple): Optional tuple (lower_limit, upper_limit) to filter outliers.
        title (str): The title of the figure.

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: The figure and axes.
    """
    fig, axes = plt.subplots(len(effect_sizes), 1, figsize=(5, 8), dpi=300)

    for i, effect_size in enumerate(effect_sizes):
        ax = axes[i]
        data = df[effect_size].dropna()

        # remove outliers beyond limits if limits provided
        if outlier_limits:
            lower_limit, upper_limit = outlier_limits
            data = data[(data > lower_limit) & (data < upper_limit)]

        ax.hist(data, bins=100, color="skyblue", edgecolor="black")
        ax.set_xlabel(effect_size, fontsize=6)
        ax.set_ylabel("Frequency", fontsize=6)
        ax.grid(ls="--", alpha=0.5)
        ax.vlines(0, *ax.get_ylim(), color="red", linestyle="--", linewidth=1)

        divider = make_axes_locatable(ax)
        box_ax = divider.append_axes("top", size="20%", pad=0.01, sharex=ax)
        box_ax.boxplot(
            data, vert=False, patch_artist=True, boxprops=dict(facecolor="lightgray")
        )
        box_ax.axis("off")
        for outlier in box_ax.findobj(match=plt.Line2D):
            outlier.set_markersize(3)
            outlier.set_alpha(0.1)

        # optional: log scale if necessary
        if max([p.get_height() for p in ax.patches]) > 10:
            ax.set_yscale("log")

        ax.tick_params(axis="both", which="major", labelsize=6)

    plt.suptitle(title, fontsize=8)
    plt.tight_layout()


def plot_effect_size_grid(
    results_df: pd.DataFrame,
    rate_types: list,
    x_var: str,
    y_vars: list[str],
    col_titles: list[str] = None,
    figure_title: str = None,
    figsize: tuple[float] = (10, 8),
    dpi: int = 300,
    s: float = 1,
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Create a grid of plots for effect sizes.

    Args:
        results_df (pandas.DataFrame): The dataframe containing the data to plot
        rate_types (list): List of rate types to plot rows for
        x_var (str): Column name to use as x variable for each grid column
        y_vars (list[str]): List of column names to use as y variables for each grid column
        col_titles (list[str], optional): Titles for each column of plots
        figure_title (str, optional): Title for the overall figure
        figsize (tuple[float], optional): Size of the figure (width, height)
        dpi (int, optional): Resolution of the figure
        s (float, optional): Size of scatter points

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
    """
    fig, axes = plt.subplots(len(rate_types), len(y_vars), figsize=figsize, dpi=dpi)

    if len(y_vars) == 1:  # if only one column, convert axes to 2D array
        axes = axes.reshape(-1, 1)

    # share x axis within columns
    for col in range(len(y_vars)):
        for i in range(1, len(rate_types)):
            axes[i, col].sharex(axes[0, col])

    for i, rate_type in enumerate(rate_types):
        rate_df = results_df[results_df["st_calcification_unit"] == rate_type]

        # Add rate type label only to first column (for the whole row)
        display_name = plot_config.RATE_TYPE_MAPPING.get(rate_type, rate_type)
        axes[i][0].set_ylabel(display_name, fontsize=6)

        for j, y_var in enumerate(y_vars):
            axes[i][j].scatter(rate_df[x_var], rate_df[y_var], s=s)

    # add column titles to the first row only
    if col_titles:
        for j, col_title in enumerate(col_titles):
            axes[0][j].set_title(col_title, fontsize=6)

    # format axes
    for ax in axes.flatten():
        ax.grid(visible=True, alpha=0.3)
        # zero effect line
        ax.axhline(y=0, color="gray", linestyle="--", alpha=0.7, linewidth=0.8)
        ax.tick_params(axis="both", which="major", labelsize=6)
        ax.ticklabel_format(
            style="scientific", axis="y", scilimits=(-2, 2), useMathText=True
        )
        ax.yaxis.get_offset_text().set_fontsize(6)

    # clear x tick labels for all but the bottom row
    for i in range(len(rate_types) - 1):
        for j in range(len(y_vars)):
            axes[i, j].tick_params(axis="x", labelbottom=False)
    for j in range(len(y_vars)):
        axes[-1, j].set_xlabel(x_var, fontsize=6)

    if figure_title:
        plt.suptitle(figure_title, fontsize=10)
    plt.tight_layout()

    return fig, axes


def create_faceted_dotplot_with_percentages(
    df: pd.DataFrame,
    top_n: int = 10,
    groupby: str = "taxa",
    omission_threshold: int = 10,
    width: float = 5,
    height: float = 10,
    dpi: int = 300,
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes] | matplotlib.figure.Figure:
    """
    Create a faceted dotplot with percentages for the top N species in each taxonomic group.

    Args:
        df (pd.DataFrame): The dataframe containing the data.
        top_n (int): The number of top species to plot.
        groupby (str): The column to group by.
        omission_threshold (int): The threshold for omitting species.
        width (float): The width of the figure.
        height (float): The height of the figure.
        dpi (int): The dpi of the figure.

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes] | matplotlib.figure.Figure: The figure and axes if requested via return_fig.
    """
    count_df = df.groupby(["species", "genus", groupby])["n"].sum().reset_index()

    # get taxa that meet the omission threshold (too few samples to be shown clearly)
    group_counts = count_df.groupby(groupby)["n"].sum()
    unique_in_group = sorted(
        [
            group
            for group in group_counts.index
            if group_counts[group] >= omission_threshold
        ]
    )
    n_in_group = len(unique_in_group)
    # order by number of samples
    unique_in_group = sorted(
        unique_in_group, key=lambda x: group_counts[x], reverse=True
    )

    fig, axes = plt.subplots(
        1, n_in_group, figsize=(width * n_in_group, height), sharey=False, dpi=dpi
    )

    if n_in_group == 1:
        axes = [axes]

    for i, group in enumerate(unique_in_group):
        ax = axes[i]

        group_data = count_df[count_df[groupby] == group].copy()
        total_count = group_data["n"].sum()

        group_data.n = group_data.n.astype(int)
        topn = group_data.nlargest(top_n, "n")
        other_species = group_data[~group_data["species"].isin(topn["species"])]

        if len(other_species) > 0:
            other_count = other_species["n"].sum()

            other_sum = pd.DataFrame(
                {
                    "species": [f"Other ({len(other_species)} species)"],
                    "genus": ["Various"],
                    groupby: [group],
                    "n": [other_count],
                }
            )

            plot_data = pd.concat([topn, other_sum], ignore_index=True)
        else:
            plot_data = topn

        plot_data["species_label"] = plot_data.apply(
            lambda row: f"{row['species']} ({(row['n'] / total_count * 100):.1f}%)",
            axis=1,
        )

        plot_data = plot_data.sort_values("n", ascending=True)

        unique_genera = sorted(plot_data["genus"].unique())
        genus_palette = dict(
            zip(unique_genera, sns.color_palette("Spectral", len(unique_genera)))
        )
        # specify 'various' genus as black
        genus_palette["Various"] = "black"

        colors = [genus_palette[genus] for genus in plot_data["genus"]]
        y_positions = range(len(plot_data))
        ax.scatter(plot_data["n"], y_positions, c=colors, s=100)

        # add count labels
        for j, (_, row) in enumerate(plot_data.iterrows()):
            ax.text(
                row["n"]
                + 0.02 * ax.get_xlim()[1]
                + len(str(row["n"])) * 0.02 * ax.get_xlim()[1],
                j,
                f"{row['n']}",
                va="center",
            )
        for j, (_, row) in enumerate(plot_data.iterrows()):
            ax.plot([0, row["n"]], [j, j], "gray", alpha=0.3)

        # formatting
        ax.set_yticks(y_positions)
        ax.set_yticklabels(plot_data["species_label"])
        max_count = plot_data["n"].max()
        ax.set_xlim(0, max_count * 1.3)  # extra space for the count labels
        ax.grid(axis="x", linestyle="--", alpha=0.7)
        ax.set_xlabel("Count", fontsize=12)
        if i == 0:
            ax.set_ylabel("Species", fontsize=12)
        # legend
        legend_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=genus_palette[genus],
                markersize=10,
            )
            for genus in unique_genera
        ]

        total_species = len(group_data)
        ax.set_title(
            f"{group.capitalize() if group != 'CCA' else 'CCA'}\n(Total: {int(total_count)} samples, {total_species} species)",
            fontsize=14,
        )
        ax.legend(legend_handles, unique_genera, title="Genus", loc="lower right")

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    return fig, axes


def plot_core_grouping_scatter(
    calcification_data_df: pd.DataFrame,
    global_anomaly_df: pd.DataFrame,
    yaxis_limits: tuple[float] = (-250, 250),
    figsize: tuple[float] = (12, 14),
    dpi: int = 150,
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Plot scatter plots of st_relative_calcification vs delta_t and delta_ph for each core_grouping, with point color and size encoded by the other variable and standard error, highlighting reasonable forecasted climatologies.

    Args:
        calcification_data_df (pd.DataFrame): Dataframe with core_grouping, delta_t, delta_ph, st_relative_calcification and variance.
        global_anomaly_df (pd.DataFrame): Climatology dataframe for threshold calculation.
        yaxis_limits (tuple): Y-axis limits for all plots.
        figsize (tuple): Figure size.
        dpi (int): Dots per inch for figure quality.

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: Matplotlib figure and axes array.
    """
    (_, max_t, min_ph, _) = climatology_processing.calculate_extreme_climatology_values(
        global_anomaly_df
    )
    core_groupings = calcification_data_df["core_grouping"].unique()
    n_groups = len(core_groupings)
    fig, axes = plt.subplots(
        n_groups,
        2,
        figsize=figsize,
        sharey=True,
        dpi=dpi,
    )

    plt.subplots_adjust(wspace=0.3)

    scatter_objects = [None, None]  # capture the scatter plots for color bars

    # order core_groupings by number of samples
    core_groupings = calcification_data_df["core_grouping"].value_counts().index
    for cg_index, cg in enumerate(core_groupings):
        cg_ax = axes[cg_index, :]
        cg_df = calcification_data_df[calcification_data_df["core_grouping"] == cg]

        # point sizes - handle zero or negative variance by clipping to min positive
        st_var = cg_df["st_relative_calcification_var"].clip(lower=1e-5)
        size = 200 * 1 / np.sqrt(st_var)

        # right column: delta_t vs st_relative_calcification, colored by delta_ph
        scatter_right = cg_ax[1].scatter(
            cg_df["delta_t"],
            cg_df["st_relative_calcification"],
            c=cg_df["delta_ph"],
            s=size,
            alpha=0.5,
            edgecolor="black",
            linewidth=0.5,
            vmin=calcification_data_df["delta_ph"].min(),
            vmax=calcification_data_df["delta_ph"].max(),
        )

        # left column: delta_ph vs st_relative_calcification, colored by delta_t
        scatter_left = cg_ax[0].scatter(
            cg_df["delta_ph"],
            cg_df["st_relative_calcification"],
            c=cg_df["delta_t"],
            s=size,
            alpha=0.5,
            edgecolor="black",
            linewidth=0.5,
            vmin=calcification_data_df["delta_t"].min(),
            vmax=calcification_data_df["delta_t"].max(),
        )

        # store the first scatter for colorbars
        if cg_index == 0:
            scatter_objects[0] = scatter_left
            scatter_objects[1] = scatter_right

        for i, a in enumerate(cg_ax.flatten()):
            a.set_ylim(yaxis_limits)
            a.vlines(
                0, yaxis_limits[0], yaxis_limits[1], color="darkgray", linestyle="--"
            )
            a.hlines(0, -1, 6, color="darkgray", linestyle="--")
            if i == 1:
                if cg_index == n_groups - 1:
                    a.set_xlabel("Temperature anomaly ($^\\circ$C)")
                a.set_xlim(-1, 6)
                climatology_plotting.shade_forecast_boundaries(a, "temp", max_t)
                a.text(
                    0.98,
                    0.98,
                    f"# of samples = {len(cg_df)}\n# of studies = {cg_df['original_doi'].nunique()}",
                    ha="right",
                    va="top",
                    transform=a.transAxes,
                    fontsize=18,
                )
            else:
                if cg_index == n_groups - 1:
                    a.set_xlabel("pH anomaly (Total scale)")
                a.text(
                    0.02,
                    0.98,
                    str(cg),
                    ha="left",
                    va="top",
                    transform=a.transAxes,
                    fontsize=18,
                )
                climatology_plotting.shade_forecast_boundaries(a, "ph", min_ph)
                a.set_xlim(-1, 0.1)

    # add horizontal color bars at the bottom of each column
    fig.colorbar(
        scatter_objects[0],
        ax=axes[:, 0],
        orientation="horizontal",
        shrink=0.8,
        pad=0.15,
        label="Temperature anomaly (°C)",
        anchor=(0.5, -1.6),
    )

    fig.colorbar(
        scatter_objects[1],
        ax=axes[:, 1],
        orientation="horizontal",
        shrink=0.8,
        pad=0.5,
        label="pH anomaly (Total scale)",
        anchor=(0.5, -1.6),
    )

    shaded_patch = mpatches.Patch(
        color="gray",
        alpha=0.4,
        label="Experimental treatments beyond reasonable forecasted climatologies",
    )

    fig.legend(
        handles=[shaded_patch],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=1,
        fontsize=16,
        frameon=False,
    )

    fig.text(
        -0.03,
        0.5,
        r"Relative calcification rate $\left(100 \times \frac{\mathrm{rate_{treatment}} - \mathrm{rate_{control}}}{\mathrm{rate_{control}}}\right)$",
        va="center",
        rotation="vertical",
        fontsize=18,
    )

    plt.tight_layout()
    return fig, axes


def plot_effect_sizes_summary(
    calcification_data_df: pd.DataFrame,
    annotate_axes_with_letters: bool = True,
    dpi: int = 300,
    figsize: tuple[float] = (8, 7),
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Combined visualization of the number of effect sizes per core grouping and number of samples.

    Args:
        calcification_data_df (pd.DataFrame): DataFrame that must contain the columns
            'core_grouping' and 'n'.

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: The figure and axes.
    """
    colours = plot_config.CG_COLOURS
    fig, axes = plt.subplots(2, 1, figsize=figsize, dpi=dpi, sharex=False)

    # Plot 1: Number of effect sizes per core grouping
    count_order = calcification_data_df["core_grouping"].value_counts().index
    sns.countplot(
        data=calcification_data_df,
        x="core_grouping",
        order=count_order,
        palette=colours,
        hue="core_grouping",
        ax=axes[0],
        width=0.5,
    )
    for p in axes[0].patches:
        axes[0].annotate(
            f"{int(p.get_height())}",
            (p.get_x() + p.get_width() / 2.0, p.get_height()),
            ha="center",
            va="bottom",
            fontsize=10,
            color="black",
            xytext=(0, 5),
            textcoords="offset points",
        )
    axes[0].set_ylim(
        0, calcification_data_df["core_grouping"].value_counts().max() * 1.1
    )
    axes[0].set_ylabel("Number of effect sizes")
    axes[0].grid(ls="--", alpha=0.5)
    axes[0].set_xlabel(" ")

    # Plot 2: Number of effect sizes per core grouping by number of samples
    cg_sample_counts = (
        calcification_data_df.groupby("core_grouping")["n"]
        .sum()
        .sort_values(ascending=False)
    )
    sns.barplot(
        x=cg_sample_counts.index,
        y=cg_sample_counts.values,
        palette=colours,
        hue=cg_sample_counts.index,
        ax=axes[1],
        width=0.5,
    )
    for p in axes[1].patches:
        axes[1].annotate(
            f"{int(p.get_height())}",
            (p.get_x() + p.get_width() / 2.0, p.get_height()),
            ha="center",
            va="bottom",
            fontsize=10,
            color="black",
            xytext=(0, 5),
            textcoords="offset points",
        )
    axes[1].set_ylim(0, cg_sample_counts.max() * 1.1)
    axes[1].set_xlabel("Core Grouping")
    axes[1].set_ylabel("Number of Samples")
    axes[1].grid(ls="--", alpha=0.5)

    # annotate subplots with letters
    if annotate_axes_with_letters:
        for i, ax in enumerate(axes.flatten()):
            ax.annotate(
                chr(65 + i),
                xy=(0.97, 0.92),
                xycoords="axes fraction",
                fontsize=18,
                ha="right",
                va="center",
            )

    plt.tight_layout()
    return fig, axes


def plot_ordered_effect_size_bar(
    clim_filtered_data_df: pd.DataFrame,
    effect_col="st_relative_calcification",
    ylim: tuple[float] = (-1000, 1000),
    dpi: int = 500,
    figsize: tuple[float] = (8, 4),
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Plot a bar chart of ordered (descending) effect sizes colored by core_grouping.

    Args:
        clim_filtered_data_df (pd.DataFrame): DataFrame with at least 'core_grouping' and effect_col columns.
        effect_col (str): The column name to plot as the bar height.
        ylim (tuple): Limits for the y-axis.
        dpi (int): Dots per inch for figure.
        figsize (tuple): Figure size.

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: The figure and axes.
    """
    # order data descendingly by effect size
    ordered_df = (
        clim_filtered_data_df.sort_values(by=effect_col, ascending=False)
        .reset_index(drop=True)
        .copy()
    )

    # plot
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    colors = ordered_df["core_grouping"].map(plot_config.CG_COLOURS)
    ax.bar(
        np.arange(len(ordered_df)),
        ordered_df[effect_col],
        color=colors,
    )
    ax.hlines(0, *ax.get_xlim(), color="black", linestyle="--", linewidth=1)
    ax.set_ylim(ylim)
    ax.set_xticks([])
    ax.set_ylabel(
        "Relative calcification rate"
        "\n"
        r"$\left(100 \times \frac{\mathrm{rate_{treatment}} - \mathrm{rate_{control}}}{\mathrm{rate_{control}}}\right)$"
    )
    ax.grid(ls="--", alpha=0.5)
    ax.set_xlim(
        -len(ordered_df.index) * 0.01,
        len(ordered_df.index) * 1.01,
    )
    ax.tick_params(axis="y", labelsize=8)

    # legend for core_grouping (excluding Foraminifera)
    handles = [
        mpatches.Patch(color=color, label=group)
        for group, color in plot_config.CG_COLOURS.items()
        if group != "Foraminifera"
    ]
    ax.legend(handles=handles, title="Core Grouping", loc="upper right")

    # vertical line at first index where effect size <= 0
    try:
        first_zero_idx = ordered_df[ordered_df[effect_col] <= 0].index[0]
        ax.vlines(
            first_zero_idx,
            *ax.get_ylim(),
            color="black",
            linestyle="--",
        )
    except IndexError:
        pass

    return fig, ax


def plot_boxplots_by_groupings(
    clim_filtered_data_df: pd.DataFrame,
    effect_type="st_relative_calcification",
    figsize: tuple[float] = (8, 7),
    dpi: int = 300,
    xlim: tuple[float] = (-1000, 1000),
    show: bool = True,
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Plots three boxplots:
    1. Overall effect size (Panel A)
    2. Distribution by core_grouping (Panel B)
    3. Distribution by treatment (Panel C)

    Args:
        clim_filtered_data_df (pd.DataFrame): Pre-processed dataframe, Foraminifera already excluded if desired.
        effect_type (str): Which effect size column to plot.
        figsize (tuple): Figure size.
        dpi (int): Dots per inch for figure.
        xlim (tuple): X-axis limits for all plots.
        show (bool): If True, call plt.show().

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: The figure and axes.
    """
    # add 'overall' category for both core_grouping and treatment
    overall_df = clim_filtered_data_df.copy()
    overall_df = overall_df.assign(core_grouping="overall", treatment="overall")

    plot_df_core = clim_filtered_data_df.copy()
    plot_df_treatment = clim_filtered_data_df.copy()
    plot_df_overall = overall_df.copy()

    fig, axes = plt.subplots(
        3,
        1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={"height_ratios": [0.2, 1, 0.6]},
        dpi=dpi,
    )

    # create palettes
    palette = {group: color for group, color in plot_config.CG_COLOURS.items()}
    palette["overall"] = "grey"

    # Panel A: "Overall"
    sns.boxplot(
        capwidths=0.2,
        flierprops={"marker": "x", "alpha": 0.3},
        data=plot_df_overall,
        x=effect_type,
        y="core_grouping",
        ax=axes[0],
        palette={"overall": "grey"},
        order=["overall"],
        hue="core_grouping",
        legend=False,
    )
    axes[0].axvline(0, color="black", linestyle="--")
    axes[0].set_ylabel("")
    axes[0].set_yticks([0])
    axes[0].set_yticklabels(["Overall"])
    axes[0].axvline(
        plot_df_overall[effect_type].mean(),
        color="black",
        label="Mean effect size",
    )

    # add annotation for number of trials and studies (Panel A)
    group_df = plot_df_overall
    n_samples = group_df["n"].sum()
    n_studies = (
        group_df["original_doi"].nunique()
        if "original_doi" in group_df.columns
        else group_df["doi"].nunique()
    )
    axes[0].text(
        0.98,
        0.5,
        f"Trials: {int(n_samples)}\nStudies: {n_studies}",
        transform=axes[0].transAxes,
        ha="right",
        va="center",
        fontsize=9,
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
    )

    # Panel B: By Core Grouping
    if hasattr(plot_config, "CG_COLOURS"):
        core_grouping_order = [
            k for k in plot_config.CG_COLOURS.keys() if k != "Foraminifera"
        ]
    else:
        core_grouping_order = sorted(clim_filtered_data_df["core_grouping"].unique())

    sns.boxplot(
        capwidths=0.2,
        flierprops={"marker": "x", "alpha": 0.3},
        data=plot_df_core,
        x=effect_type,
        y="core_grouping",
        ax=axes[1],
        palette=palette,
        order=core_grouping_order,
        hue="core_grouping",
        legend=False,
    )
    axes[1].axvline(0, color="black", linestyle="--")
    axes[1].set_ylabel("")
    axes[1].axvline(
        plot_df_core[effect_type].mean(),
        color="black",
        label="Mean effect size",
    )

    # add annotation for each core_grouping
    for i, core_group in enumerate(core_grouping_order):
        group_df = plot_df_core[plot_df_core["core_grouping"] == core_group]
        n_samples = group_df["n"].sum()
        n_studies = (
            group_df["original_doi"].nunique()
            if "original_doi" in group_df.columns
            else group_df["doi"].nunique()
        )
        y_pos = i
        axes[1].text(
            0.98,
            1 - ((y_pos + 0.5) / len(core_grouping_order)),  # top to bottom
            f"Samples: {int(n_samples)}\nStudies: {n_studies}",
            transform=axes[1].transAxes,
            ha="right",
            va="center",
            fontsize=9,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
        )

    # Panel C: By Treatment
    treatment_order = ["phtot", "temp", "temp_phtot"]
    treatment_order = [
        t for t in treatment_order if t in clim_filtered_data_df["treatment"].unique()
    ]
    other_treatments = [
        t
        for t in clim_filtered_data_df["treatment"].unique()
        if t not in treatment_order
    ]
    treatment_order += other_treatments

    treatment_palette = {
        t: plot_config.TREATMENT_COLOURS.get(t, "Set3") for t in treatment_order
    }

    sns.boxplot(
        capwidths=0.2,
        flierprops={"marker": "x", "alpha": 0.3},
        data=plot_df_treatment,
        x=effect_type,
        y="treatment",
        ax=axes[2],
        palette=treatment_palette,
        order=treatment_order,
        hue="treatment",
        legend=False,
    )

    # rewrite y-axis treatment labels
    axes[2].set_yticks(list(range(len(treatment_order))))
    axes[2].set_yticklabels(
        [plot_config.TREATMENT_MAPPING.get(t, t) for t in treatment_order]
    )
    axes[2].axvline(0, color="black", linestyle="--")
    axes[2].set_xlabel(
        r"Relative calcification rate $\left(100 \times \frac{\mathrm{rate_{treatment}} - \mathrm{rate_{control}}}{\mathrm{rate_{control}}}\right)$"
    )
    axes[2].set_ylabel("")
    axes[2].axvline(
        plot_df_treatment[effect_type].mean(),
        color="black",
        label="Mean effect size",
    )

    # add annotation for each treatment
    for i, treatment in enumerate(treatment_order):
        group_df = plot_df_treatment[plot_df_treatment["treatment"] == treatment]
        n_samples = group_df["n"].sum()
        n_studies = (
            group_df["original_doi"].nunique()
            if "original_doi" in group_df.columns
            else group_df["doi"].nunique()
        )
        y_pos = i
        axes[2].text(
            0.98,
            1 - ((y_pos + 0.5) / len(treatment_order)),
            f"Trials: {int(n_samples)}\nStudies: {n_studies}",
            transform=axes[2].transAxes,
            ha="right",
            va="center",
            fontsize=9,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
        )

    for ax in axes.flatten():
        ax.xaxis.grid(alpha=0.5, ls="--")
        ax.set_xlim(*xlim)
    plot_utils.annotate_axes_with_letters(axes, xy=(25, 15))

    fig.legend(
        [
            Line2D([0], [0], color="black", linestyle="--"),
            Line2D([0], [0], color="black", linestyle="-"),
        ],
        ["Zero effect size", "Mean effect size"],
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncols=2,
    )

    plt.tight_layout()
    if show:
        plt.show()
    return fig, axes
