# general
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.markers import MarkerStyle
from shapely.geometry import MultiPolygon, Polygon

from calcification_meta_analysis.plotting import plot_utils


def site_distribution_spatial_plot(
    calcification_df: pd.DataFrame,
    bioerosion_df: pd.DataFrame,
    all_reefs: gpd.GeoDataFrame,
    tropical_realms: gpd.GeoDataFrame,
    legend_order: list[str],
    grid_size: float = 2,
    min_marker_size: float = 20,
    max_marker_size: float = 200,
    figsize: tuple[float, float] = (15, 15),
    dpi: int = 300,
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Master spatial plot function for study density and reef/realm overlays.

    Args:
        calcification_df (pd.DataFrame): DataFrame of calcification study locations.
        bioerosion_df (pd.DataFrame): DataFrame of bioerosion study locations.
        all_reefs (gpd.GeoDataFrame): Reef polygons for the background.
        tropical_realms (gpd.GeoDataFrame): Polygon/multipolygon geometries for ecoregional realms.
        legend_order (list): List of realm names (in order) for legend and color assignment.
        grid_size (float): Grid and marker scaling parameters.
        min_marker_size (float): Grid and marker scaling parameters.
        max_marker_size (float): Grid and marker scaling parameters.
        figsize (tuple): Figure size.
        dpi (int): DPI.

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: The figure and axes.
    """

    def plot_study_density(
        df: pd.DataFrame,
        color: str,
        label: str,
        ax: matplotlib.axes.Axes,
        grid_size: int = 2,
        min_marker_size: int = 20,
        max_marker_size: int = 200,
        return_grid: bool = False,
    ) -> tuple[pd.DataFrame, int, int, int, int] | None:
        """
        Plot study density with improved overlapping marker visualization.

        Args:
            df (pd.DataFrame): The dataframe containing the study locations.
            color (str): The color of the study type.
            label (str): The label of the study type.
            ax (matplotlib.axes.Axes): The axes to plot the study density on.
            grid_size (int): The size of the grid cells.
            min_marker_size (int): The minimum marker size.
            max_marker_size (int): The maximum marker size.
            return_grid (bool): Whether to return the grid data.

        Returns:
            tuple[pd.DataFrame, int, int, int, int] | None: The grid data and the minimum and maximum counts and marker sizes.
                If return_grid is True, the grid data is returned, otherwise None is returned.
        """
        df = df.copy()
        # ensure coordinates are numeric and drop nans
        df["latitude"] = pd.to_numeric(df["latitude"], errors="coerce")
        df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")
        df = df[pd.notna(df["latitude"]) & pd.notna(df["longitude"])]
        # bin coordinates into grid cells specified by grid_size argument
        lat_bins = (df["latitude"] // grid_size) * grid_size
        lon_bins = (df["longitude"] // grid_size) * grid_size
        grid = (
            df.groupby([lat_bins, lon_bins])["original_doi"]
            .nunique()
            .reset_index(name="count")
        )
        grid.columns = ["lat_bin", "lon_bin", "count"]
        # scale marker size based on number of studies in each grid cell
        if not grid["count"].empty:
            min_count, max_count = grid["count"].min(), grid["count"].max()
            if min_count == max_count:
                sizes = [max_marker_size for _ in grid["count"]]
            else:
                sizes = min_marker_size + (grid["count"] - min_count) / (
                    max_count - min_count
                ) * (max_marker_size - min_marker_size)
        else:
            sizes = []

        # plot split markers for overlapping circles indicating bioerosion and calcification studies
        ax.scatter(
            grid["lon_bin"] + grid_size / 2,
            grid["lat_bin"] + grid_size / 2,
            marker=MarkerStyle("o", fillstyle="right" if color == "blue" else "left"),
            s=sizes,
            color=color,
            alpha=0.7,
            linewidth=0,
            transform=ccrs.PlateCarree(),
            zorder=105,  # back half of each point (always z=105)
        )

        ax.scatter(
            grid["lon_bin"] + grid_size / 2,
            grid["lat_bin"] + grid_size / 2,
            marker=MarkerStyle("o", fillstyle="left" if color == "blue" else "right"),
            s=sizes,
            color=color,
            alpha=0.7,
            linewidth=0,
            transform=ccrs.PlateCarree(),
            zorder=110,  #  front half of each point (always z=110)
        )
        ax.scatter(
            grid["lon_bin"] + grid_size / 2,
            grid["lat_bin"] + grid_size / 2,
            s=sizes,
            facecolors="none",
            edgecolors="white",
            alpha=0.7,
            zorder=120,
            linewidth=0.5,
            transform=ccrs.PlateCarree(),
        )  # plot empty circle with white outline

        if return_grid:
            return grid, min_count, max_count, min_marker_size, max_marker_size
        else:
            return None

    def get_size_for_count(
        count: int,
        min_count: int,
        max_count: int,
        min_marker_size: int,
        max_marker_size: int,
    ) -> float:
        """Get the size for a given count based on the min and max count and marker size.

        Args:
            count (int): The count to get the size for.
            min_count (int): The minimum count.
            max_count (int): The maximum count.
            min_marker_size (int): The minimum marker size.
            max_marker_size (int): The maximum marker size.

        Returns:
            float: The size for the given count.
        """
        if min_count == max_count:
            return max_marker_size
        else:
            return min_marker_size + (count - min_count) / (max_count - min_count) * (
                max_marker_size - min_marker_size
            )

    # prepare colormap
    cmap = plt.get_cmap("viridis")
    norm = mcolors.Normalize(vmin=0, vmax=len(legend_order) - 1)
    fig, ax = plt.subplots(
        subplot_kw={"projection": ccrs.PlateCarree()}, dpi=dpi, figsize=figsize
    )
    # plot reefs in the background
    all_reefs.plot(ax=ax, color="sandybrown", alpha=0.5, linewidth=0, zorder=-5)

    # plot outlines oftropical realms
    for idx, row in tropical_realms.iterrows():
        geom = row.geometry
        if isinstance(geom, Polygon):
            exteriors = [geom.exterior]
        elif isinstance(geom, MultiPolygon):
            exteriors = [poly.exterior for poly in geom.geoms]
        else:
            continue
        # assign color based on legend order to get nice gradient
        color_idx = legend_order.index(idx) if idx in legend_order else 0
        color = plt.get_cmap("viridis")(color_idx / (len(legend_order) - 1))
        for ext in exteriors:
            xs, ys = ext.xy
            ax.plot(xs, ys, color=color, linewidth=2, alpha=0.5)

    # plot study densities
    _, bio_min_count, bio_max_count, bio_min_marker_size, bio_max_marker_size = (
        plot_study_density(
            bioerosion_df,
            color="blue",
            label="Bioerosion",
            ax=ax,
            grid_size=grid_size,
            min_marker_size=min_marker_size,
            max_marker_size=max_marker_size,
            return_grid=True,
        )
    )
    (
        _,
        calc_min_count,
        calc_max_count,
        calc_min_marker_size,
        calc_max_marker_size,
    ) = plot_study_density(
        calcification_df,
        color="red",
        label="Calcification",
        ax=ax,
        grid_size=grid_size,
        min_marker_size=min_marker_size,
        max_marker_size=max_marker_size,
        return_grid=True,
    )

    # add legend for the study types (ignoring marker size) using dummy points
    calc_dummy = plt.Line2D(
        [],
        [],
        marker="o",
        color="w",
        markerfacecolor="red",
        markersize=8,
        label="Calcification",
    )
    bio_dummy = plt.Line2D(
        [],
        [],
        marker="o",
        color="w",
        markerfacecolor="blue",
        markersize=8,
        label="Bioerosion",
    )

    legend1 = ax.legend(
        handles=[calc_dummy, bio_dummy],
        loc="upper left",
        fontsize=8,
        title="Study type",
    )
    ax.add_artist(legend1)

    # legend for marker sizes (study count)
    all_min_count = min(calc_min_count, bio_min_count)
    all_max_count = max(calc_max_count, bio_max_count)
    all_min_marker_size = min(calc_min_marker_size, bio_min_marker_size)
    all_max_marker_size = max(calc_max_marker_size, bio_max_marker_size)

    # representative counts
    if all_max_count - all_min_count > 2:
        legend_counts = [
            all_min_count,
            (all_min_count + all_max_count) // 2,
            all_max_count,
        ]
        legend_counts = sorted(set(int(x) for x in legend_counts))
    else:
        legend_counts = sorted(set(int(x) for x in [all_min_count, all_max_count]))

    # create legend handles
    size_handles = [
        plt.scatter(
            [],
            [],
            s=get_size_for_count(
                c,
                all_min_count,
                all_max_count,
                all_min_marker_size,
                all_max_marker_size,
            ),
            color="gray",
            alpha=0.7,
            edgecolor="white",
            linewidth=0.5,
        )
        for c in legend_counts
    ]
    size_labels = [f"{c} study" if c == 1 else f"{c} studies" for c in legend_counts]

    legend2 = ax.legend(
        size_handles,
        size_labels,
        scatterpoints=1,
        loc="upper right",
        fontsize=8,
        title="Density",
        frameon=True,
        borderpad=0.8,
        labelspacing=0.8,
        handletextpad=1.2,
    )
    ax.add_artist(legend2)

    # prepare legend handles for realm outlines
    realm_patches = []
    for i, realm in enumerate(legend_order):
        color = cmap(norm(i))
        patch = mpatches.Patch(
            facecolor="none",
            edgecolor=color,
            linewidth=2,
            label=realm,
        )
        realm_patches.append(patch)

    # add legend for realm outlines beneath plot
    ax.legend(
        realm_patches,
        legend_order,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=len(legend_order),
    )

    ax = plot_utils.format_geo_axes(ax)
    return fig, ax
