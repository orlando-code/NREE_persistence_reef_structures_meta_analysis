# general
import unicodedata
from pathlib import Path

import cbsyst.helpers as cbh
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from calcification_meta_analysis.utils import config, file_ops, utils


### helper functions
def cluster_values(values: list, tolerance: float) -> list:
    """
    Cluster values based on their proximity.

    Args:
        values (array-like): Values to cluster.
        tolerance (float): Tolerance for clustering.

    Returns:
        list: List of clusters, where each cluster is a list of values.
    """
    if len(values) == 0:  # return empty list if no values
        return []

    # sort values
    sorted_values = np.sort(values)

    # initialize first cluster
    clusters = [[sorted_values[0]]]

    # cluster remaining values
    for val in sorted_values[1:]:
        # check if value is sufficiently close to the last value in the current cluster
        if np.abs(val - np.mean(clusters[-1])) < tolerance:
            # add to current (most recent) cluster
            clusters[-1].append(val)
        else:  # if not close enough
            # start new cluster
            clusters.append([val])

    return clusters


def aggregate_treatments_with_individual_samples(df: pd.DataFrame) -> pd.DataFrame:
    """For treatments with only one sample (most often those for which raw, sample-level data was extracted), aggregate the data to get means and standard deviations of the treatment groups."""
    aggregated_df = (
        df.groupby(
            [
                "doi",
                "species_types",
                "treatment_level_ph",
                "treatment_level_t",
                "calcification_unit",
            ]
        )
        .filter(lambda group: (group["n"] == 1).all())
        .groupby(
            [
                "doi",
                "species_types",
                "treatment_level_ph",
                "treatment_level_t",
                "calcification_unit",
            ]
        )
        .agg(
            ecoregion=("ecoregion", "first"),  # metadata
            lat_zone=("lat_zone", "first"),
            latitude=("latitude", "first"),
            longitude=("longitude", "first"),
            location=("location", "first"),
            taxa=("taxa", "first"),
            genus=("genus", "first"),
            species=("species", "first"),
            family=("family", "first"),
            core_grouping=("core_grouping", "first"),
            authors=("authors", "first"),
            year=("year", "first"),
            treatment_group=("treatment_group", "first"),
            treatment=("treatment", "first"),
            dic=("dic", "mean"),  # carbonate chemistry
            dic_sd=("dic", "std"),
            pco2=("pco2", "mean"),
            pco2_sd=("pco2", "std"),
            phtot=("phtot", "mean"),
            phtot_sd=("phtot", "std"),
            temp=("temp", "mean"),
            temp_sd=("temp", "std"),
            sal=("sal", "mean"),
            sal_sd=("sal", "std"),
            irr=("irr", "mean"),  # irradiance
            irr_sd=("irr", "std"),
            calcification=("calcification", "mean"),  # calcification
            calcification_sd=("calcification", "std"),
            st_calcification=("st_calcification", "mean"),
            st_calcification_sd=("st_calcification", "std"),
            st_calcification_unit=("st_calcification_unit", "first"),
            n=("n", "count"),
        )
        .reset_index()
    )
    # remove rows with n=1
    df_no_ones = df[df["n"] != 1]
    # append the aggregated data to the DataFrame
    df_no_ones = pd.concat([df_no_ones, aggregated_df], ignore_index=True)
    return df_no_ones


def calc_sd_from_se(se: float, n: int) -> float:
    """Calculate standard deviation from standard error and sample size

    Args:
        se (float): standard error
        n (int): number of samples

    Returns:
        float: standard deviation
    """
    return se * np.sqrt(n)


### raw file wrangling
def preprocess_df(df: pd.DataFrame) -> pd.DataFrame:
    """Clean dataframe fields and standardise for future processing"""
    ### basic cleaning
    df.columns = df.columns.str.normalize("NFKC").str.replace(
        "μ", "u"
    )  # replace any unicode versions of 'μ' with 'u'
    df = df.map(
        lambda x: unicodedata.normalize("NFKD", str(x)).replace("\xa0", " ")
        if isinstance(x, str)
        else x
    )  # clean non-breaking spaces from string cells
    # general processing
    df.rename(
        columns=file_ops.read_yaml(config.resources_dir / "mapping.yaml")[
            "sheet_column_map"
        ],
        inplace=True,
    )  # rename columns to agree with cbsyst output
    df.columns = (
        df.columns.str.lower()
    )  # columns lower case headers for less confusing access later on
    df.columns = df.columns.str.replace(
        " ", "_"
    )  # process columns to replace whitespace with underscore
    df.columns = df.columns.str.replace(
        "[()]", "", regex=True
    )  # remove '(' and ')' from column names
    df["year"] = pd.to_datetime(
        df["year"], format="%Y"
    )  # datetime format for later plotting

    ### deal with duplicate dois: flag up duplicate dois which also have 'include' as 'yes'
    inclusion_df = df[df["include"] == "yes"]
    duplicate_dois = inclusion_df[inclusion_df.duplicated(subset="doi", keep=False)]
    if not duplicate_dois.empty and not all(pd.isna(duplicate_dois["doi"])):
        print("\nDuplicate DOIs found, treat with caution:")
        print([doi for doi in duplicate_dois.doi.unique() if doi is not np.nan])

    ### formating: fill down necessary repeated metadata values
    df[["doi", "year", "authors", "location", "species_types", "taxa"]] = (
        df[["doi", "year", "authors", "location", "species_types", "taxa"]]
        .infer_objects(copy=False)
        .ffill()
    )
    df[["coords", "cleaned_coords"]] = df.groupby("doi")[
        ["coords", "cleaned_coords"]
    ].ffill()  # fill only as far as the next DOI

    ### deal with missing n
    if (
        df["n"].dtype == "object"
    ):  # Only perform string operations if column contains strings
        df = df[
            ~df["n"].str.contains("~", na=False)
        ]  # remove any rows in which 'n' has '~' in the string
        df = df[df.n != "M"]  # remove any rows in which 'n' is 'M'

    ### infer data types
    df.loc[:, df.columns != "year"] = df.loc[:, df.columns != "year"].apply(
        utils.safe_to_numeric
    )
    problem_cols = [
        "irr",
        "ipar",
        "sal",
    ]  # some columns have rogue strings when they should all contain numbers: in this case, convert unconvertable values to NaN
    for col in problem_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # remove any columns with 'unnamed' in the header: these are an artefact from messing around outside the spreadsheets necessary columns
    df = df.loc[:, ~df.columns.str.contains("^unnamed")]
    return df


### carbonate chemistry
def calculate_carb_chem(row: pd.Series, out_values: list) -> dict:
    """(Re)calculate carbonate chemistry parameters from the dataframe row and return a dictionary."""
    try:
        out_dict = cbh.Csys(
            pHtot=row["phtot"],
            TA=row["ta"],
            T_in=row["temp"],
            S_in=35 if pd.isna(row["sal"]) else row["sal"],
        )
        out_dict = {
            key.lower(): value for key, value in out_dict.items()
        }  # lower the keys of the dictionary to ensure case-insensitivity
        return {
            key: (
                out_dict.get(key.lower(), None)[0]
                if isinstance(out_dict.get(key.lower()), (list, np.ndarray))
                else out_dict.get(key.lower(), None)
            )
            for key in out_values
        }
    except Exception as e:
        print(f"Error: {e}")


### assigning treatment groups
def determine_control_conditions(df: pd.DataFrame) -> dict:
    """Identify the rows corresponding to min temperature and/or max pH.

    Args:
        df (pd.DataFrame): Input dataframe with columns 'doi', 'temp', 'phtot', etc.

    Returns:
        dict: Dictionary with control conditions for each treatment group.
    """
    grouped = df.groupby("treatment_group")

    control_treatments = {}

    for group, sub_df in grouped:
        group = int(group)  # convert group to integer for semantics
        min_temp = (
            sub_df.loc[sub_df["temp"].idxmin()]["temp"]
            if not any(sub_df["phtot"].isna())
            else None
        )  # Row with minimum temperature
        max_pH = (
            sub_df.loc[sub_df["phtot"].idxmax()]["phtot"]
            if not any(sub_df["phtot"].isna())
            else None
        )  # Row with maximum pH

        control_treatments[group] = {
            "control_t_in": min_temp,
            "control_phtot": max_pH,
        }

    return control_treatments


def assign_treatment_groups(
    df: pd.DataFrame,
    control_T: float,
    control_pH: float,
    t_mapping: dict,
    ph_mapping: dict,
    irr_group: float,
) -> pd.DataFrame:
    """Assign treatment groups based on temperature and pH values.

    Args:
        df (pd.DataFrame): Input dataframe with columns 'doi', 'temp', 'phtot', etc.
        control_T (float): Control temperature value.
        control_pH (float): Control pH value.
        t_mapping (dict): Mapping of temperature values to cluster indices.
        ph_mapping (dict): Mapping of pH values to cluster indices.
        irr_group (float): Irradiance group identifier.

    Returns:
        pd.DataFrame: Dataframe with added treatment group columns.
    """
    # apply classification to each row in this group
    for idx in df.index:
        row = df.loc[idx]

        # get temperature cluster level (0 is control)
        t_level = None
        if not np.isnan(row["temp"]) and control_T is not None:
            t_cluster_idx = t_mapping.get(row["temp"])
            control_cluster_idx = t_mapping.get(control_T)
            if t_cluster_idx is not None and control_cluster_idx is not None:
                t_level = t_cluster_idx - control_cluster_idx

        # get pH cluster level (0 is control)
        ph_level = None
        if not np.isnan(row["phtot"]) and control_pH is not None:
            ph_cluster_idx = ph_mapping.get(row["phtot"])
            control_cluster_idx = ph_mapping.get(control_pH)
            if ph_cluster_idx is not None and control_cluster_idx is not None:
                ph_level = (
                    control_cluster_idx - ph_cluster_idx
                )  # reverse order since higher pH is control

        # determine clusters for cases where there is only one of t or ph
        if t_level is None and ph_level is not None:
            t_level = 0
        if ph_level is None and t_level is not None:
            ph_level = 0

        # determine if values are in control clusters   # TODO: not currently capturing rare case when studies have both T and pH varied from control with no intermediary values
        is_control_T = t_level == 0 if t_level is not None else False
        is_control_pH = ph_level == 0 if ph_level is not None else False

        # classify the treatment
        if is_control_T and is_control_pH:
            treatment = "cTcP"
        elif is_control_T:
            treatment = "cTtP"
        elif is_control_pH:
            treatment = "tTcP"
        elif not (is_control_T or is_control_pH):
            treatment = "tTtP"
        else:
            treatment = "uncertain"

        # Update the treatment info in the result dataframe
        df.loc[idx, "treatment_group"] = treatment
        df.loc[idx, "treatment_level_t"] = t_level if t_level is not None else np.nan
        df.loc[idx, "treatment_level_ph"] = ph_level if ph_level is not None else np.nan
        df.loc[idx, "irr_group"] = irr_group

    return df


### climatology
def process_climatology_csv(fp: str, index_col: str = "doi") -> pd.DataFrame:
    df = pd.read_csv(fp).drop(columns=["data_ID", "Unnamed: 0"])
    # rename columns to be less wordy
    df = (
        (df.copy())
        .replace({"2021_2040": 2030, "2041_2060": 2050, "2081_2100": 2090})
        .infer_objects(copy=False)
    )

    return df.set_index(index_col) if index_col else df


def convert_climatology_csv_to_multiindex(
    fp: str, locations_yaml_fp: str
) -> pd.DataFrame:
    """
    Convert the climatology CSV file to a multi-index DataFrame.
    """
    df = process_climatology_csv(fp, index_col="doi")  # load the CSV file

    var = "ph" if "ph" in str(fp.name) else "sst" if "sst" in str(fp.name) else None
    if not var:
        raise ValueError(
            "File path must contain 'ph' or 'sst' to determine variable type."
        )
    df = pd.concat(
        [
            df.iloc[:, :4],
            df.iloc[:, 4:].rename(
                columns=lambda col: col if var in col else f"{var}_{col}"
            ),
        ],
        axis=1,
    )

    # load locations yaml as dataframe
    locations_df = pd.DataFrame(file_ops.read_yaml(locations_yaml_fp)).T
    # reorder columns to be latitude, longitude, location
    locations_df = locations_df[["latitude", "longitude", "location"]]

    # merge locations with sst_df
    df = df.merge(
        locations_df,
        left_index=True,
        right_index=True,
        how="left",
        suffixes=("", "_right"),
    )
    # remove any columns ending with '_right' or the 'index_right' column specifically
    columns_to_drop = [
        col for col in df.columns if col.endswith("_right") or col == "index_right"
    ]
    df = df.drop(columns=columns_to_drop)
    df.reset_index(inplace=True, names="doi")

    return df


def generate_location_specific_anomalies(df: pd.DataFrame, scenario_var: str = "sst"):
    """Generate location-specific anomalies for each location in the dataframe"""
    df = (
        df.sort_index()
    )  # sort the index to avoid PerformanceWarning about lexsort depth
    locations = df.index.unique()
    anomaly_rows = []  # to hold newmods inputs
    metadata_rows = []  # to track what each row corresponds to

    for location in tqdm(
        locations, desc=f"Generating batched anomalies for {scenario_var}"
    ):
        location_df = df.loc[location]
        scenarios = location_df["scenario"].unique()

        for scenario in scenarios:
            scenario_df = location_df[location_df["scenario"] == scenario]
            time_frames = [1995] + list(scenario_df.time_frame.unique())

            for time_frame in time_frames:
                if time_frame == 1995:
                    base = scenario_df[
                        f"mean_historical_{scenario_var}_30y_ensemble"
                    ].mean()
                    mean_scenario = base - base
                    p10_scenario = (
                        scenario_df[
                            f"percentile_10_historical_{scenario_var}_30y_ensemble"
                        ].mean()
                        - base
                    )
                    p90_scenario = (
                        scenario_df[
                            f"percentile_90_historical_{scenario_var}_30y_ensemble"
                        ].mean()
                        - base
                    )
                else:
                    time_scenario_df = scenario_df[
                        scenario_df["time_frame"] == time_frame
                    ]
                    mean_scenario = time_scenario_df[
                        f"mean_{scenario_var}_20y_anomaly_ensemble"
                    ].mean()
                    p10_scenario = time_scenario_df[
                        f"{scenario_var}_percentile_10_anomaly_ensemble"
                    ].mean()
                    p90_scenario = time_scenario_df[
                        f"{scenario_var}_percentile_90_anomaly_ensemble"
                    ].mean()
                    # Generate predictions for mean, p10, and p90 scenarios
                for percentile, anomaly in [
                    ("mean", mean_scenario),
                    ("p10", p10_scenario),
                    ("p90", p90_scenario),
                ]:
                    anomaly_rows.append([anomaly])
                    metadata_rows.append(
                        {
                            "doi": location[0],
                            "location": location[1],
                            "longitude": location[2],
                            "latitude": location[3],
                            "scenario_var": scenario_var,
                            "scenario": scenario,
                            "time_frame": time_frame,
                            # 'anomaly_value': anomaly,
                            "percentile": percentile,
                        }
                    )
    return pd.concat(
        [
            pd.DataFrame(metadata_rows),
            pd.DataFrame(anomaly_rows, columns=["anomaly_value"]),
        ],
        axis=1,
    )


def process_emissions_sheet(sheet_df: pd.DataFrame, scenario_name: str) -> pd.DataFrame:
    """Process the emissions sheet DataFrame as provided from supplementary data of https://doi.org/10.5194/gmd-13-3571-2020."""
    sheet_df = sheet_df[["Gas", "CO2"]].iloc[3:]  # years labelled 'Gas'
    sheet_df.rename(columns={"Gas": "year", "CO2": scenario_name}, inplace=True)
    sheet_df["year"] = pd.to_numeric(sheet_df["year"], errors="coerce")
    return sheet_df


def load_and_format_coral_cover(
    filepath: Path = config.data_dir / "coral_cover_locs.xlsx",
) -> pd.DataFrame:
    """Load coral cover data and reshape for analysis."""
    df = (
        pd.read_excel(filepath, sheet_name="Sheet1")
        .infer_objects()
        .astype("float64", errors="raise")
        .drop_duplicates()
    )
    # melt and extract features
    id_vars = ["Reef", "ReefLat", "ReefLon"]
    df_long = df.melt(id_vars=id_vars, var_name="metric", value_name="value")
    df_long[["year", "coral_type", "scenario"]] = df_long["metric"].str.extract(
        r"(\d{4})\s+(Mounding|Branching)\s+(rcp\d{2})"
    )
    df_long = df_long.drop(columns="metric")
    df_long["year"] = df_long["year"].astype(int)
    # standardize column names
    df_long = df_long.rename(
        columns={"ReefLat": "lat", "ReefLon": "lon", "Reef": "reef_id"}
    )
    df_long.columns = df_long.columns.str.lower()
    df_long = df_long.apply(
        lambda col: col.str.lower() if col.dtype == "object" else col
    )
    # wide table for comparisons between predicted coral cover in2050/2100
    wide_df = df_long.pivot_table(
        index=["reef_id", "scenario", "coral_type", "lat", "lon"],
        columns="year",
        values="value",
        aggfunc="first",
    )
    return wide_df


def compute_cover_difference(wide_df: pd.DataFrame) -> pd.DataFrame:
    """Compute cover change between 2050 and 2100 and order by mounding coral loss."""
    df_2050 = wide_df[2050].reset_index()
    df_2100 = wide_df[2100].reset_index()
    df_diff = pd.merge(
        df_2050,
        df_2100,
        on=["reef_id", "scenario", "lat", "lon", "coral_type"],
        suffixes=("_2050", "_2100"),
    )
    df_diff["diff"] = df_diff[2100] - df_diff[2050]
    # map labels of RCP to SSP for consistency (this is fairly approximate since rcp and ssp are not exactly the same)
    rcp_ssp_map = {"rcp85": "ssp585", "rcp45": "ssp245", "rcp26": "ssp126"}
    df_diff["scenario"] = df_diff["scenario"].map(rcp_ssp_map)
    # add scenario-wide std, merge in
    diff_std = df_diff.groupby("scenario")["diff"].std().reset_index()
    df_diff = pd.merge(df_diff, diff_std, on="scenario", suffixes=("", "_std"))
    # sort reefs by mounding mean
    mounding_means = (
        df_diff[df_diff["coral_type"] == "mounding"]
        .groupby("reef_id", as_index=False)["diff"]
        .mean()
    )
    sorted_reef_ids = mounding_means.sort_values(by="diff", ascending=False)["reef_id"]
    # adjust for plotting
    df_diff["diff"] = df_diff["diff"] * 100
    df_diff["reef_id_cat"] = pd.Categorical(
        df_diff["reef_id"], categories=sorted_reef_ids, ordered=True
    )
    df_diff["cover_2100"] = df_diff[2100] * 100
    return df_diff
