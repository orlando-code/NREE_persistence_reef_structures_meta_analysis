# general
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# stats
import statsmodels.api as sm
from matplotlib.legend_handler import HandlerBase
from scipy.stats import norm
from tqdm.auto import tqdm

from calcification_meta_analysis.analysis import analysis_utils, metafor
from calcification_meta_analysis.plotting import plot_config
from calcification_meta_analysis.processing import climatology as climatology_processing


# --- core analysis calculations ---
def calc_relative_rate(
    mu1: float,
    mu2: float,
    sd1: float = None,
    sd2: float = None,
    n1: int = None,
    n2: int = None,
    epsilon: float = 1e-6,
) -> tuple[float, float]:
    """
    Calculate percent change between two means with error propagation.

    Examples:
    - 10 to 20: +100%
    - 10 to 0: -100%
    - 10 to -10: -200%

    Args:
        mu1, mu2 (float):   Mean values to compare (mu1=reference/baseline, mu2=new value)
        se1, se2 (float, optional):   Standard errors of mu1 and mu2
        epsilon (float, optional):  Small value to stabilize calculations when means are close to zero

    Returns:
        pc (float): Percent change
        se_pc (float or None):  Standard error of the percent change
    """
    se1, se2 = sd1 / np.sqrt(n1), sd2 / np.sqrt(n2)
    if mu1 == 0 and mu2 == 0:  # special case: both means are exactly zero
        pc = 0  # no change between means

        if (
            se1 is not None and se2 is not None
        ):  # if SEs are provided, calculate uncertainty
            # when both means are zero, consider the ratio of SEs to estimate uncertainty
            # this represents how much percentage change we would expect if the values
            # fluctuated by ±1 SE from zero
            if (
                se1 > 0
            ):  # scale based on the potential percentage fluctuations around zero
                se_pc = 100 * se2 / se1
            else:  # if se1 is zero but se2 is not, technically infinite uncertainty
                se_pc = float("inf") if se2 > 0 else 0
            return pc, se_pc
        return pc

    if mu1 == 0:  # special case: baseline is zero
        raise ValueError("Baseline is zero")

    # standard percent change calculation
    pc = ((mu2 - mu1) / abs(mu1)) * 100

    # return only PC if no standard errors provided
    if se1 is None or se2 is None:
        return pc

    # error propagation via partial derivatives
    dpc_dmu1 = (-mu2 / (mu1**2)) * 100
    dpc_dmu2 = (1 / abs(mu1)) * 100
    var_pc = (dpc_dmu1**2 * se1**2) + (dpc_dmu2**2 * se2**2)

    return pc, var_pc


def calc_absolute_rate(
    mu1: float,
    mu2: float,
    sd1: float = None,
    sd2: float = None,
    n1: int = None,
    n2: int = None,
) -> tuple[float, float]:
    """Calculate the simple difference between two means with error propagation.

    Args:
        mu1, mu2 (float):   Mean values to compare (mu1=reference/baseline, mu2=new value)
        sd1, sd2 (float, optional):   Standard deviations of mu1 and mu2
        n1, n2 (int): number of samples in group 1 (control) and group 2 (treatment)

    Returns:
        tuple: absolute difference between means, standard error of the difference
    """
    abs_diff = mu2 - mu1

    # if standard deviations are provided, calculate the uncertainty
    if sd1 is not None and sd2 is not None and n1 is not None and n2 is not None:
        se1 = sd1 / np.sqrt(n1)
        se2 = sd2 / np.sqrt(n2)

        # error propagation - calculate partial derivatives
        d_abs_diff_dmu1 = -1
        d_abs_diff_dmu2 = 1

        # calculate standard error using error propagation
        var_abs_diff = (d_abs_diff_dmu1**2 * se1**2) + (d_abs_diff_dmu2**2 * se2**2)

        return abs_diff, var_abs_diff

    return abs_diff, None


def calc_bias_correction(n1: int, n2: int) -> float:
    """Calculate bias correction for Cohen's d metric: https://www.campbellcollaboration.org/calculator/equations

    Args:
        n1, n2 (int): number of samples in group 1 (control) and group 2 (treatment)

    Returns:
        float: bias correction factor
    """
    return 1 - 3 / (4 * (n1 + n2 - 2) - 1)


def calc_cohens_d(
    mu1: float, mu2: float, sd1: float, sd2: float, n1: int, n2: int
) -> tuple[float, float]:
    """Calculate Cohen's d metric: https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/hedgeg.htm

    Args:
        mu1, mu2 (float):   Mean values to compare (mu1=reference/baseline, mu2=new value)
        sd1, sd2 (float, optional):   Standard deviations of mu1 and mu2
        n1, n2 (int): number of samples in group 1 and group 2

    Returns:
        tuple[float, float]: Cohen's d and its variance
    """
    sd_pooled = calc_pooled_sd(n1, n2, sd1, sd2)
    d = (mu2 - mu1) / sd_pooled if sd_pooled != 0 else 0
    d_var = (n1 + n2) / (n1 * n2) + d**2 / (2 * (n1 + n2))
    return d, d_var


def calc_pooled_sd(n1: int, n2: int, sd1: float, sd2: float) -> float:
    """Calculate pooled standard deviation for two groups.
    N.B. BH (2021) uses simple average

    Args:
        n1, n2 (int): number of samples in group 1 (control) and group 2 (treatment)
        sd1, sd2 (float, optional):   Standard deviations of mu1 and mu2

    Returns:
        float: pooled standard deviation
    """
    return np.sqrt(((n1 - 1) * sd1**2 + (n2 - 1) * sd2**2) / (n1 + n2 - 2))


def calc_hedges_g(
    mu1: float, mu2: float, sd1: float, sd2: float, n1: int, n2: int
) -> tuple[float, float]:
    """Calculate Hedges G metric: https://www.campbellcollaboration.org/calculator/equations

    Args:
        mu1, mu2 (float):   Mean values to compare (mu1=reference/baseline, mu2=new value)
        sd1, sd2 (float, optional):   Standard deviations of mu1 and mu2
        n1, n2 (int): number of samples in group 1 and group 2

    Returns:
        float: Hedges G metric
    """
    d, d_var = calc_cohens_d(mu1, mu2, sd1, sd2, n1, n2)
    bias_correction = calc_bias_correction(n1, n2)

    hg = d * bias_correction
    hg_var = d_var * bias_correction**2
    return hg, hg_var


# --- meta-analysis functions ---
def calc_cooks_distance(data: pd.Series) -> pd.Series:
    """
    Calculate Cook's distance for a given data series.
    """
    # if data is not numeric
    if not pd.api.types.is_numeric_dtype(data):
        # convert data to numeric
        data = pd.to_numeric(data, errors="coerce")

    # fit OLS model
    X = sm.add_constant(np.asarray(data))
    try:
        model = sm.OLS(data, X).fit()
    except ValueError:
        # convert data to numeric if it is not already
        model = sm.OLS(data, X).fit()
    # calculate Cook's distance
    influence = model.get_influence()
    cooks_d = influence.cooks_distance[0]

    return cooks_d


def calc_cooks_threshold(
    data: pd.Series, nparams: int, threshold_type: str = "liberal"
) -> float:
    """
    Calculate the Cook's distance threshold for outlier detection via a reproducible numerical replacement of eyeballing for outliers in the distance-study graph.

    Args:
        data (pd.Series): The data to calculate the threshold for.
        nparams (int): The number of parameters in the model.
        threshold_type (str): The type of threshold to calculate.

    Returns:
        float: The Cook's distance threshold.
    """
    n = len(data)
    print(n, nparams)
    if threshold_type == "conservative":
        threshold = 4 / (n - nparams)
    elif threshold_type == "liberal":
        threshold = 2 * np.sqrt(nparams / (n - nparams - 1))
    else:
        raise ValueError(f"Invalid threshold type: {threshold_type}")
    return threshold


def remove_cooks_outliers(
    df: pd.DataFrame,
    effect_type: str = "hedges_g",
    nparams: int = 3,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Remove outliers from a DataFrame based on Cook's distance.
    """
    data = df.copy()
    # calculate cooks distance
    cooks_threshold = calc_cooks_threshold(data[effect_type], nparams=nparams)
    # calculate cooks distance
    data["cooks_d"] = calc_cooks_distance(data[effect_type])

    # remove outliers
    data_no_outliers = data[data["cooks_d"] < cooks_threshold]
    outliers = data[data["cooks_d"] >= cooks_threshold]
    print(
        f"\nRemoved {len(outliers)} outlier(s) (from {len(data)} samples) based on Cook's distance threshold of {cooks_threshold:.2f}"
    ) if verbose else None
    return data_no_outliers, outliers


# --- calculating treatment effects ---
def calculate_effect_for_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate Hedges' g for a DataFrame of experimental results.

    Args:
        df: DataFrame containing experimental data

    Returns:
        pandas.DataFrame: DataFrame with calculated effect sizes
    """
    # copy to avoid modifying original
    result_df = df.copy()

    effect_cols = [
        "delta_t",
        "delta_ph",
        "cohens_d",
        "cohens_d_var",
        "hedges_g",
        "hedges_g_var",
        "relative_calcification",
        "relative_calcification_var",
        "absolute_calcification",
        "absolute_calcification_var",
        "st_relative_calcification",
        "st_relative_calcification_var",
        "st_absolute_calcification",
        "st_absolute_calcification_var",
    ]
    for col in effect_cols:
        result_df[col] = np.nan

    # remove any rows with n=1 (since will have no error)
    n_1_rows = result_df[result_df["n"] == 1].shape[0]
    print(f"Removing {n_1_rows} {'row' if n_1_rows == 1 else 'rows'} with n=1")
    result_df = result_df[result_df["n"] != 1]  # this shortens the df by 1

    # group by relevant factors and apply processing
    grouped_data = []
    doi_bar = tqdm(result_df.doi.unique())
    for doi in doi_bar:
        doi_bar.set_description(f"Calculating effect sizes for {doi}")
        study_df = result_df[result_df["doi"] == doi]
        for _, irr_df in study_df.groupby("irr_group"):
            for _, species_df in irr_df.groupby("species_types"):
                df = calculate_row_effect(species_df)
                if isinstance(df, pd.Series):
                    df = pd.DataFrame([df].T)
                if df is not None:
                    grouped_data.append(df)
    if isinstance(grouped_data, list):
        valid_dfs = [
            df
            for df in grouped_data
            if df is not None and not df.empty and not df.isna().all().all()
        ]
        if valid_dfs:
            df = (
                pd.concat(valid_dfs).sort_index().copy()
            )  # sort index to avoid lexsort depth warning and create a copy to avoid fragmentation
        else:
            # return empty DataFrame with same columns and dtypes as expected output
            df = pd.DataFrame(
                columns=df.columns
                if len(grouped_data) > 0 and grouped_data[0] is not None
                else None
            )
    df = df.sort_values(by="doi").copy().reset_index(drop=True)
    df["ID"] = df.index

    df.loc[df.phtot.isna(), "delta_ph"] = (
        0  # assign delta_ph = 0 where phtot is NaN (assumes this variable was controlled throughout the experiment)
    )
    df.loc[df.temp.isna(), "delta_t"] = (
        0  # similarly, assign delta_t = 0 where temp is NaN
    )

    # replace any 0 values in "*_var" columns with mean of that column (there's no such thing as no error)
    for col in df.columns:
        if col.endswith("_var"):
            mean_value = df[col].mean()
            df[col] = df[col].replace(0, mean_value).infer_objects(copy=False)

    return df


def calculate_control_values(control_df: pd.DataFrame) -> pd.Series:
    """Calculates the representative control series by averaging numeric columns."""
    if control_df.empty:
        return pd.Series(dtype=object)  # return empty series if no control data
    if len(control_df) > 1:
        numeric_cols = control_df.select_dtypes(include="number").columns
        # create a Series with first values for non-numeric, means for numeric
        control_series = control_df.iloc[0].copy()
        for col in numeric_cols:
            if col == "n":
                control_series[col] = control_df[col].sum(skipna=True)
            else:
                control_series[col] = control_df[col].mean(skipna=True)
    else:
        control_series = control_df.iloc[0].copy()
    return control_series


def calculate_row_effect(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the effect size for each row of a dataframe which contains a control and some number of treatments.
    This method ignores treatment groups, other than to identify the control.

    Args:
        df (pd.DataFrame): DataFrame containing control and treatment data

    Returns:
        pd.DataFrame: Row with calculated effect sizes and additional metadata
    """
    # identify control
    control_df = df[df["treatment"] == "control"]
    if control_df.empty:
        print(
            f"Control dataframe is empty for DOI: {df.doi.iloc[0]}. Consider assigning manually (currently losing {df.shape[0]} rows)"
        )
        return None
    control_series = calculate_control_values(control_df)
    # calculate effect size for each row in treatment_df and create a list of results
    effect_rows = []
    for _, row in df.iterrows():
        if row["treatment"] == "control":
            continue
        effect_row = calc_treatment_effect_for_row(row, control_series)
        effect_rows.append(effect_row)
    if effect_rows:
        return pd.concat(effect_rows, axis=1).T.copy()
    return None


def aggregate_by_treatment_group(df: pd.DataFrame) -> pd.Series:
    """
    Aggregate a DataFrame by treatment group. Useful for when samples are individual datapoints, or multiple slightly-different controls are sent.

    Args:
        df: DataFrame containing data for a specific treatment group

    Returns:
        pandas.Series: Series containing aggregated data
    """
    aggregation = df.agg({"calcification": ["mean", "std"], "n": "sum"})
    control_row = df.iloc[0].copy()
    control_row["calcification"] = aggregation["calcification"]["mean"]
    control_row["calcification_sd"] = aggregation["calcification"]["std"]
    control_row["n"] = aggregation["n"]["count"]
    return control_row


def calc_treatment_effect_for_row(
    treatment_row: pd.Series, control_data: pd.Series
) -> pd.Series:
    """
    Calculate the effect size (Hedges' g or relative calcification) and append additional columns for a treatment row.

    Args:
        treatment_row: Row containing treatment data
        control_data: Dictionary containing control group data

    Returns:
        pandas.Series: Row with calculated effect sizes and additional metadata
    """
    # raw values
    mu_t, sd_t, n_t = (
        treatment_row["calcification"],
        treatment_row["calcification_sd"],
        treatment_row["n"],
    )
    mu_c, sd_c, n_c = (
        control_data["calcification"],
        control_data["calcification_sd"],
        control_data["n"],
    )
    # standardised values
    s_mu_t, s_sd_t, _ = (
        treatment_row["st_calcification"],
        treatment_row["st_calcification_sd"],
        treatment_row["n"],
    )
    s_mu_c, s_sd_c, _ = (
        control_data["st_calcification"],
        control_data["st_calcification_sd"],
        control_data["n"],
    )
    t_in_c, ph_c = control_data["temp"], control_data["phtot"]

    if (
        np.isnan(mu_t)
        or np.isnan(mu_c)
        or np.isnan(sd_t)
        or np.isnan(sd_c)
        or np.isnan(s_mu_t)
        or np.isnan(s_mu_c)
        or np.isnan(s_sd_t)
        or np.isnan(s_sd_c)
    ):
        print(
            f"Missing data for effect size calculation. "
            f"n_t: {n_t:.3f}, n_c: {n_c:.3f} at \n[index {treatment_row.name} DOI {treatment_row['doi']}]"
            f"Raw: mu_t: {mu_t:.3f}, mu_c: {mu_c:.3f}, sd_t: {sd_t:.3f}, sd_c: {sd_c:.3f}, "
            f"Std: s_mu_t: {s_mu_t:.3f}, s_mu_c: {s_mu_c:.3f}, s_sd_t: {s_sd_t:.3f}, s_sd_c: {s_sd_c:.3f}, "
        )

    row_copy = treatment_row.copy()  # create a copy to avoid SettingWithCopyWarning

    d_effect, d_var = calc_cohens_d(mu_c, mu_t, sd_c, sd_t, n_c, n_t)  # Cohen's d
    hg_effect, hg_var = calc_hedges_g(mu_c, mu_t, sd_c, sd_t, n_c, n_t)  # Hedges' g

    rc_effect, rc_var = calc_relative_rate(mu_c, mu_t, sd_c, sd_t, n_c, n_t)

    abs_effect, abs_var = calc_absolute_rate(
        mu_c, mu_t, sd_c, sd_t, n_c, n_t
    )  # absolute differences

    st_d_effect, st_d_var = calc_cohens_d(
        s_mu_c, s_mu_t, s_sd_c, s_sd_t, n_c, n_t
    )  # standardised cohen's d
    st_hg_effect, st_hg_var = calc_hedges_g(
        s_mu_c, s_mu_t, s_sd_c, s_sd_t, n_c, n_t
    )  # standardised hedges' g

    # relative differences between standardised calcification rates
    st_rc_effect, st_rc_var = calc_relative_rate(
        s_mu_c, s_mu_t, s_sd_c, s_sd_t, n_c, n_t
    )
    # absolute differences between standardised calcification rates
    st_abs_effect, st_abs_var = calc_absolute_rate(
        s_mu_c, s_mu_t, s_sd_c, s_sd_t, n_c, n_t
    )

    # assign effect sizes
    row_copy.update(
        {
            "cohens_d": d_effect,
            "cohens_d_var": d_var,
            "hedges_g": hg_effect,
            "hedges_g_var": hg_var,
            "relative_calcification": rc_effect,
            "relative_calcification_var": rc_var,
            "absolute_calcification": abs_effect,
            "absolute_calcification_var": abs_var,
            "st_relative_calcification": st_rc_effect,
            "st_relative_calcification_var": st_rc_var,
            "st_absolute_calcification": st_abs_effect,
            "st_absolute_calcification_var": st_abs_var,
            "st_cohens_d": st_d_effect,
            "st_cohens_d_var": st_d_var,
            "st_hedges_g": st_hg_effect,
            "st_hedges_g_var": st_hg_var,
        }
    )

    # calculate metadata
    row_copy["control_temp"] = control_data["temp"]
    row_copy["treatment_temp"] = treatment_row["temp"]
    row_copy["delta_t"] = row_copy["temp"] - t_in_c
    row_copy["control_phtot"] = control_data["phtot"]
    row_copy["treatment_phtot"] = treatment_row["phtot"]
    row_copy["delta_ph"] = row_copy["phtot"] - ph_c
    row_copy["treatment_val"] = (
        row_copy["temp"] if row_copy["treatment"] == "temp" else row_copy["phtot"]
    )
    row_copy["control_calcification"] = mu_c
    row_copy["control_calcification_sd"] = sd_c
    row_copy["treatment_calcification"] = mu_t
    row_copy["treatment_calcification_sd"] = sd_t
    row_copy["st_control_calcification"] = s_mu_c
    row_copy["st_control_calcification_sd"] = s_sd_c
    row_copy["st_treatment_calcification"] = s_mu_t
    row_copy["st_treatment_calcification_sd"] = s_sd_t
    row_copy["treatment_n"] = n_t
    row_copy["control_n"] = n_c

    return row_copy


def weighted_mean_ci(
    x: np.ndarray, meas_se: np.ndarray, mu0: float = 0.0, alpha: float = 0.05
) -> tuple[float, float, float, float]:
    """Calculate the weighted mean, confidence interval, and significance level from data using inverse of variance as weights.

    Args:
        x: array-like of data
        meas_se: array-like of measurement standard errors
        mu0: null hypothesis value
        alpha: significance level

    Returns:
        tuple: weighted mean, lower confidence interval, upper confidence interval, p-value
    """
    # calculate mean weighted by inverse of variance
    x = np.asarray(x, float)
    meas_se = np.asarray(meas_se, float)
    var = meas_se**2
    w = 1.0 / var
    mu_hat = (w * x).sum() / w.sum()
    var_mu = 1.0 / w.sum()
    se_mu = np.sqrt(var_mu)
    # calculate critical z-value for significance
    z_crit = norm.ppf(1 - alpha / 2)
    z = (mu_hat - mu0) / se_mu
    # calculate p-value
    p = 2 * (1 - norm.cdf(abs(z)))
    return mu_hat, mu_hat - z_crit * se_mu, mu_hat + z_crit * se_mu, p


def summarise_group(
    df, group_col, n_col="n", doi_col="doi", effect_type="st_relative_calcification"
):
    """Summarise a group of data by calculating the mean, confidence interval, and significance."""
    rows = []
    for group, subset in df.groupby(group_col):
        data = subset[effect_type].values.astype(float)
        meas_se = subset[effect_type + "_var"].values.astype(float)
        mean, low, high, p = weighted_mean_ci(data, meas_se)

        n_trials = subset[n_col].sum() if n_col in subset.columns else len(subset)
        n_studies = (
            subset["original_doi"].nunique()
            if "original_doi" in subset.columns
            else subset[doi_col].nunique()
            if doi_col in subset.columns
            else np.nan
        )

        rows.append(
            {
                "group": group,
                "mean": mean,
                "low": low,
                "high": high,
                "n": len(subset),
                "n_trials": n_trials,
                "n_studies": n_studies,
                "p_val": p,
            }
        )
    return pd.DataFrame(rows)


def label_study_data(axis: plt.axis, row: pd.Series, x_value: float = -200) -> None:
    """Label the study data on the axis.

    Args:
        axis: The axis to label the study data on.
        row: The row of the dataframe containing the study data.
        x_value: The x-value to label the study data on.
    """
    axis.text(
        x_value,
        row["group"],
        f"Trials: {int(row['n_trials'])}\nStudies: {int(row['n_studies'])}",
        va="center",
        ha="left",
        fontsize=8,
        color="black",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
    )


class HandlerStars(HandlerBase):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        import matplotlib.text as mtext

        label = orig_handle.get_label()
        # cast to float if safe
        try:
            label = float(label)
        except ValueError:
            pass
        # invert mapping
        label = {v: k for k, v in plot_config.SIGNIFICANCE_MAPPING.items()}[label]
        # Center the text in the legend box
        t = mtext.Text(
            xdescent + width / 2,
            ydescent + height / 2,
            label,
            ha="center",
            va="center_baseline",
            fontsize=fontsize + 2,
        )
        return [t]


def construct_predicted_response_surfaces(
    df: pd.DataFrame,
    global_anomaly_df: pd.DataFrame,
    effect_type: str,
    model_formula: str = "delta_t + delta_ph - 1",
    surface_resolution: int = 1000,
    cooks_distance_type: bool = None,
):
    """
    For each core_grouping in the data, fit model and generate predicted response surfaces sampled over observed/future climatologies (global_anomaly_df).

    Args:
        df (pd.DataFrame): DataFrame of data to fit models and generate response surfaces
        global_anomaly_df (pd.DataFrame): DataFrame of global climatology anomalies
        effect_type (str): Effect type to fit models with
        surface_resolution (int): Resolution of the response surface

    Returns:
        response_df: DataFrame of response predictions for each scenario/row in global_anomaly_df (tiled for each core_grouping)
        heterogeneity_df: DataFrame with model heterogeneity stats for each core_grouping
        models_dict: Dict of fitted MetaforModel objects for each core_grouping
    """
    (_, max_t, min_ph, _) = climatology_processing.calculate_extreme_climatology_values(
        global_anomaly_df
    )

    response_df = pd.DataFrame()
    heterogeneity_df = pd.DataFrame()
    models_dict = {}

    for cg in tqdm(df["core_grouping"].unique()):
        cg_df = df[df["core_grouping"] == cg].copy()
        cg_model = metafor.MetaforModel(
            cg_df,
            effect_type=effect_type,
            formula=f"{effect_type} ~ {model_formula}",
            verbose=False,
            cooks_distance_type=cooks_distance_type,
        ).fit_model()

        models_dict[cg] = cg_model  # for downstream debugging/exploration

        # prepare model heterogeneity info
        cg_het_df = cg_model.get_heterogeneity_dataframe()
        cg_het_df["core_grouping"] = cg
        heterogeneity_df = pd.concat([heterogeneity_df, cg_het_df], axis=0)

        # predict response surface for this grouping
        model_surface, meshgrids = cg_model.predict_nd_surface_from_model(
            moderator_names=["delta_t", "delta_ph"],
            moderator_values=[
                np.linspace(0, max_t, surface_resolution),
                np.linspace(min_ph, 0, surface_resolution),
            ],
        )

        cg_climatology = global_anomaly_df.copy()
        cg_climatology.loc[:, "core_grouping"] = cg

        for k in model_surface.keys():  # for each climate scenario
            cg_climatology.loc[:, k] = cg_climatology.apply(
                lambda row: analysis_utils.get_value_from_surface(
                    model_surface[k],
                    meshgrids,
                    row["anomaly_value_sst"],
                    row["anomaly_value_ph"],
                ),
                axis=1,
            )

        # assign p-scores and certainty
        dof = len(cg_df) - len(cg_model.coefficients)
        cg_climatology = analysis_utils.assign_p_score_and_certainty(
            cg_climatology, dof=dof
        )
        response_df = pd.concat([response_df, cg_climatology], axis=0)

    return response_df, heterogeneity_df, models_dict
