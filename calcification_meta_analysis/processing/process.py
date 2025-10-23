import logging
from typing import Optional

import pandas as pd

from calcification_meta_analysis.analysis import analysis
from calcification_meta_analysis.processing import (
    carbonate_processing,
    cleaning,
    climatology,
    groups_processing,
    locations,
    processing,
)
from calcification_meta_analysis.utils import config, file_ops

logging.basicConfig(level=logging.INFO)


def process_extracted_calcification_data(
    fp: str, sheet_name: str = "all_data", selection_dict: Optional[dict] = None
) -> pd.DataFrame:
    """
    Full pipeline for processing calcification data from raw Excel to effect sizes.

    Args:
        fp (str): Path to Excel file.
        sheet_name (str): Sheet name in Excel file.
        selection_dict (dict): Optional dict for row selection.

    Returns:
        pd.DataFrame: DataFrame with effect sizes and all processing applied.
    """
    if selection_dict is None:  # exclude selected rows
        selection_dict = {"include": "yes"}

    # populate carbonate chemistry
    carbonate_df = carbonate_processing.populate_carbonate_chemistry(
        fp, sheet_name=sheet_name, selection_dict=selection_dict
    )
    treatment_group_df = process_carbonate_df_to_treatment_groups(carbonate_df)
    effect_sizes_df = process_treatment_group_df_to_effect_sizes(treatment_group_df)
    return effect_sizes_df


def process_carbonate_df_to_treatment_groups(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Process carbonate dataframe to treatment groups."""
    # assign treatment groups
    treatment_group_df = processing.assign_treatment_groups_multilevel(df)
    # drop rows with nan in treatment (this is due to no treatment conditions identified)
    treatment_group_df = treatment_group_df.dropna(subset=["treatment"])
    # aggregate treatments with individual samples
    return groups_processing.aggregate_treatments_rows_with_individual_samples(
        treatment_group_df
    )


def process_treatment_group_df_to_effect_sizes(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Process treatment group dataframe to effect sizes."""
    effect_sizes_df = analysis.calculate_effect_for_df(df)
    # infer dtypes for columns that are not numeric
    effect_sizes_df = effect_sizes_df.infer_objects()
    # calculate the dcalcification_dvariable values
    effect_sizes_df = processing.calculate_dvar(effect_sizes_df)
    return effect_sizes_df


def process_extracted_df_to_effect_sizes(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Process clean extracted dataframe to effect sizes."""
    # do the processing on the cleaned dataframe here to make it match the old version
    df = df.rename(
        columns=file_ops.read_yaml(config.resources_dir / "mapping.yaml")[
            "sheet_column_map"
        ]
    )
    df.columns = df.columns.str.lower()
    # create species_types column
    df["species_types"] = df[["genus", "species"]].agg(" ".join, axis=1)
    df["original_doi"] = df["doi"].str.split("-LOC").str[0]
    # create original_doi column
    treatment_group_df = process_carbonate_df_to_treatment_groups(df)
    # standard calcification rate
    treatment_group_df = cleaning.standardise_calcification_rates(treatment_group_df)
    effect_sizes_df = process_treatment_group_df_to_effect_sizes(treatment_group_df)

    # save locations to yaml file
    locations.save_locations_information(
        effect_sizes_df
    )  # N.B. all coords already populated in pre-processing

    # clean up the dataframe to be clean again
    # undo column name mapping
    effect_sizes_df = effect_sizes_df.rename(
        columns=file_ops.read_yaml(config.resources_dir / "mapping.yaml")[
            "sheet_column_map"
        ],
    )
    return effect_sizes_df


def process_carbonate_df_to_effect_sizes(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Process carbonate dataframe to effect sizes."""
    treatment_group_df = process_carbonate_df_to_treatment_groups(df)
    effect_sizes_df = process_treatment_group_df_to_effect_sizes(treatment_group_df)
    return effect_sizes_df


def process_climatology_data(
    experimental_df: pd.DataFrame,
    ph_clim_path: Optional[str] = None,
    sst_clim_path: Optional[str] = None,
    locations_path: Optional[str] = None,
    experiment_type: Optional[str] = "calcification",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merge processed data with climatology and compute global average anomalies.

    Args:
        experimental_df (pd.DataFrame): Processed calcification DataFrame (with latitudes and longitudes).
        ph_clim_path (str): Path to pH climatology CSV.
        sst_clim_path (str): Path to SST climatology CSV.
        locations_path (str): Path to locations YAML.

    Returns:
        pd.DataFrame: DataFrame with local anomalies by scenario and time_frame (2030, 2050, 2090).
        pd.DataFrame: DataFrame with global average anomalies by scenario and time_frame.
    """
    ph_clim_path = ph_clim_path or (
        config.climatology_data_dir / "ph_scenarios_output_table_site_locations.csv"
    )
    sst_clim_path = sst_clim_path or (
        config.climatology_data_dir / "sst_scenarios_output_table_site_locations.csv"
    )
    locations_path = (
        locations_path or (config.resources_dir / "locations.yaml")
    )  # N.B. locations.yaml has far fewer locs than all_locations.yaml. Need to check what's happened here.

    logging.info("Loading climatology data...")
    # if experiment_type == "calcification":
    ph_climatology = climatology.convert_climatology_csv_to_multiindex(
        ph_clim_path, locations_path
    )
    sst_climatology = climatology.convert_climatology_csv_to_multiindex(
        sst_clim_path, locations_path
    )

    sst_ph_climatology_df = pd.merge(sst_climatology, ph_climatology)
    # elif experiment_type == "bioerosion":
    #     sst_ph_climatology_df = pd.read_csv(
    #         config.climatology_data_dir / "site_locations_with_MMM_and_pH.csv"
    #     )

    sst_ph_climatology_df_mi = sst_ph_climatology_df.set_index(
        ["doi", "location", "longitude", "latitude"]
    )
    experimental_df_mi = experimental_df.set_index(
        ["doi", "location", "longitude", "latitude"]
    )
    local_climatology_df = experimental_df_mi.join(
        sst_ph_climatology_df_mi, how="inner"
    )

    # logging.info(
    #     f"Unique locations in climatology: {len(sst_ph_climatology_df_mi.index.unique())}, "
    #     f"locations in working experimental dataframe: {len(experimental_df.drop_duplicates('doi', keep='first'))}"
    # )

    # exclude aquaria locations # TODO: make this more robust/automated
    # local_climatology_df = local_climatology_df[
    #     ~local_climatology_df.index.get_level_values("location").str.contains(
    #         "monaco|portugal|uk", case=False, na=False
    #     )
    # ]

    ph_anomalies = climatology.generate_location_specific_climatology_anomalies(
        local_climatology_df, "ph"
    )
    sst_anomalies = climatology.generate_location_specific_climatology_anomalies(
        local_climatology_df, "sst"
    )
    # Merge the two DataFrames on the relevant columns
    merged_anomalies = pd.merge(
        ph_anomalies,
        sst_anomalies,
        on=[
            "doi",
            "location",
            "longitude",
            "latitude",
            "scenario",
            "time_frame",
            "percentile",
        ],
        suffixes=("_ph", "_sst"),
    ).drop(columns=["scenario_var_ph", "scenario_var_sst"])

    global_anomaly_df = (
        merged_anomalies.groupby(["scenario", "time_frame", "percentile"])[
            ["anomaly_value_ph", "anomaly_value_sst"]
        ]
        .mean()
        .reset_index()
    )  # average spatially
    # calculate global average anomalies for dataframe
    global_future_anomaly_df = (
        local_climatology_df.reset_index()
        .groupby(["scenario", "time_frame"])
        .agg(
            mean_sst_anomaly=("mean_sst_20y_anomaly_ensemble", "mean"),
            mean_ph_anomaly=("mean_ph_20y_anomaly_ensemble", "mean"),
        )
        .reset_index()
    )
    # return local_climatology_df, future_climatology_df
    return local_climatology_df, global_future_anomaly_df, global_anomaly_df
