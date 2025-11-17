import re

import numpy as np
import pandas as pd
import rpy2.robjects as ro
import statsmodels.api as sm
from scipy.interpolate import make_interp_spline

from calcification_meta_analysis.analysis import analysis
from calcification_meta_analysis.utils import config, file_ops


def preprocess_df_for_meta_model(
    df: pd.DataFrame,
    effect_type: str = "st_relative_calcification",
    effect_type_var: str = None,
    treatment: list[str] | str = None,
    formula_components: dict = None,
    verbose: bool = True,
    apply_approx_cooks_threshold: bool = False,
) -> pd.DataFrame:
    data = df.copy()
    df["doi"] = df["doi"].astype(str)

    # select only rows relevant to treatment
    if treatment:
        if isinstance(treatment, list):
            data = data[data["treatment"].astype(str).isin(treatment)]
        else:
            data = data[data["treatment"] == treatment]

    n_investigation = len(data)
    # remove nans for subset effect_type
    required_columns = _get_required_columns(
        effect_type, formula_components, effect_type_var
    )
    data = data.dropna(subset=required_columns)
    data = data.convert_dtypes()
    n_nans = n_investigation - len(data)

    # remove outliers
    nparams = get_number_of_params(formula_components)
    data, cooks_outliers = (
        analysis.remove_cooks_outliers(
            data, effect_type=effect_type, nparams=nparams, verbose=False
        )
        if apply_approx_cooks_threshold
        else (data, [])
    )

    if verbose:
        # summarise processing
        print("\n----- PROCESSING SUMMARY -----")
        print("Treatment: ", treatment)
        print("Total samples in input data: ", len(df))
        print("Total samples of relevant investigation: ", n_investigation)
        print("Total samples dropped due to OLS Cook's distance: ", len(cooks_outliers))
        print("Dropped due to NaN values: ", n_nans)
        nan_counts = df[required_columns].isna().sum()
        for col, count in nan_counts.items():
            if count > 0:
                print(f"\t{col}: {count} NaN(s)")
        print("Dropped due to Cook's distance: ", len(cooks_outliers))
        print(
            f"Final sample count: {len(data)} ({len(cooks_outliers) + n_nans + (len(df) - n_investigation)} rows dropped)\n"
        )

    return data.infer_objects()


def get_number_of_params(formula_components: dict) -> int:
    """Get the number of parameters in a formula."""
    return (
        len(formula_components["raw_predictors"]) + 1
        if formula_components["intercept"]
        else len(formula_components["raw_predictors"])
    )


def generate_metaregression_formula(
    effect_type: str,
    treatment: str = None,
    include_intercept: bool = False,
) -> str:
    treatment_vars = _get_treatment_vars(treatment)
    variable_mapping = file_ops.read_yaml(config.resources_dir / "mapping.yaml")[
        "meta_model_factor_variables"
    ]
    factor_vars = [f"factor({v})" for v in treatment_vars if v in variable_mapping]
    other_vars = [v for v in treatment_vars if v not in variable_mapping]

    # combine into formula string
    formula = f"{effect_type} ~ {' + '.join(other_vars + factor_vars)}"

    # remove intercept if specified
    return formula + " - 1" if not include_intercept else formula


def _get_required_columns(
    effect_type: str,
    formula_components: dict,
    effect_type_var: str = None,
    required_columns: list[str] | None = None,
):
    if not formula_components:
        raise ValueError("Formula components are required")
    if formula_components["raw_predictors"] == ["0"]:
        raise ValueError(
            "Metafor formula requires at least one predictor e.g. an intercept"
        )
    # get required columns from formula components
    formula_requirements = formula_components["raw_predictors"]

    effect_type_var = (
        f"{effect_type}_var" if effect_type_var is None else effect_type_var
    )
    base_columns = [
        "doi",
        "ID",
        "core_grouping",
        "st_calcification_unit",
        "st_control_calcification",
        "st_treatment_calcification",
        effect_type,
        effect_type_var,
    ]
    use_columns = (
        base_columns + formula_requirements
        if formula_requirements != ["1"]
        else base_columns
    )
    return use_columns + required_columns if required_columns else use_columns


def _get_treatment_vars(treatment: str) -> list[str]:
    """Get the treatment variables for the model."""

    def _get_treatment_var_from_single_treatment(treatment: str) -> list[str]:
        if treatment == "phtot":
            return ["delta_ph"]
        elif treatment == "temp":
            return ["delta_t"]
        elif treatment in ["phtot_mv", "temp_mv", "phtot_temp_mv"]:
            return ["delta_ph", "delta_t"]
        else:
            raise ValueError(f"Unknown treatment: {treatment}")

    treatment_vars = []
    if isinstance(treatment, list):
        for t in treatment:
            treatment_vars.extend(_get_treatment_var_from_single_treatment(t))
    else:
        treatment_vars.extend(_get_treatment_var_from_single_treatment(treatment))

    return list(set(treatment_vars))


def get_formula_components(formula: str) -> dict:
    """
    Extracts the response variable, predictors, and intercept flag from a formula string.
    Handles complex interactions like I(temp^2):I(phtot^2) and mixed term types.

    Args:
        formula (str): A formula string, e.g. "y ~ delta_t:delta_ph + factor(core_grouping) + I(temp^2):I(phtot^2)"

    Returns:
        dict: {
            "response": str,
            "raw_predictors": list[str],  # All unique variables used
            "linear_terms": list[str],    # Simple linear terms
            "factor_terms": list[str],    # factor() terms
            "nonlinear_terms": list[str], # I() terms
            "interaction_terms": list[str], # interaction terms (a:b)
            "factorial_terms": list[str],   # full factorial terms (a*b)
            "intercept": bool
        }
    """
    # split formula into response and predictors
    response_part, predictor_part = [p.strip() for p in formula.split("~", 1)]

    # handle intercept removal
    predictor_part = predictor_part.replace(" ", "")
    has_intercept = "-1" not in predictor_part
    predictor_part = predictor_part.replace("-1", "")

    # split on '+' to get individual terms
    raw_terms = [term.strip() for term in predictor_part.split("+") if term.strip()]

    # initialize collections for different term types
    linear_terms = []
    factor_terms = []
    nonlinear_terms = []
    interaction_terms = []
    factorial_terms = []
    raw_predictors = set()

    for term in raw_terms:
        if not term:
            continue

        # check for factorial terms (a*b expands to a + b + a:b)
        if "*" in term:
            factorial_terms.append(term)
            # extract variables from factorial term
            factors = term.split("*")
            for factor in factors:
                var = _extract_variable_name(factor.strip())
                if var:
                    raw_predictors.add(var)

        # check for interaction terms (a:b)
        elif ":" in term:
            interaction_terms.append(term)
            # extract variables from interaction
            interactors = term.split(":")
            for interactor in interactors:
                var = _extract_variable_name(interactor.strip())
                if var:
                    raw_predictors.add(var)

        # check for factor terms
        elif term.startswith("factor("):
            # extract variable name from factor(variable)
            match = re.search(r"factor\(([^)]+)\)", term)
            if match:
                raw_predictors.add(match.group(1))
                factor_terms.append(match.group(1))

        # check for nonlinear terms I(...)
        elif term.startswith("I("):
            nonlinear_terms.append(term)
            # extract variable name from I(variable^power)
            var = _extract_variable_name(term)
            if var:
                raw_predictors.add(var)

        # simple linear term
        else:
            linear_terms.append(term)
            raw_predictors.add(term)

    return {
        "response": response_part,
        "raw_predictors": sorted(list(raw_predictors)),
        "factor_terms": factor_terms,
        "linear_terms": linear_terms,
        "nonlinear_terms": nonlinear_terms,
        "interaction_terms": interaction_terms,
        "factorial_terms": factorial_terms,
        "intercept": has_intercept,
        "formula": formula,
    }


def _extract_variable_name(term: str) -> str:
    """
    Extract the core variable name from various term formats.

    Examples:
        delta_t -> delta_t
        I(delta_t^2) -> delta_t
        factor(core_grouping) -> core_grouping
    """

    # handle I(...) terms
    if term.startswith("I("):
        # extract variable from I(variable^power) or I(variable)
        match = re.search(r"I\(([^)^]+)", term)
        if match:
            return match.group(1)

    # handle factor(...) terms
    elif term.startswith("factor("):
        match = re.search(r"factor\(([^)]+)\)", term)
        if match:
            return match.group(1)

    # simple variable name
    else:
        return term.strip()

    return ""


def p_score(
    prediction: float,
    se: float,
    dof: int,
    null_value: float = 0,
) -> float:
    """
    Calculate the p-value for a given prediction and standard error.

    Args:
        prediction (float): The prediction.
        se (float): The standard error.
        dof (int): The degrees of freedom for the t-distribution (n-k where n is the number of observations and k is the number of coefficients).
        null_value (float): The null value.

    Returns:
        float: The p-value.
    """
    if se == 0:
        return 0
    t_score = (prediction - null_value) / se
    from scipy.stats import t as scipy_t

    return 2 * (
        1 - scipy_t.cdf(abs(t_score), df=dof)
    )  # two-tailed p-value with t-distribution


def assign_certainty(p_score: float) -> int:
    """
    Assign certainty levels based on p-value.
    """
    if p_score < 0.01:
        return 3  # very high certainty
    elif p_score < 0.05:
        return 2  # high certainty
    elif p_score < 0.1:
        return 1  # medium certainty
    else:
        return 0  # low certainty


def assign_p_score_and_certainty(predictions: pd.DataFrame, dof: int):
    """
    Assign p-scores and certainty to a dataframe of predictions.

    Args:
        predictions (pd.DataFrame): The dataframe containing the predictions.
        dof (int): The degrees of freedom for the t-distribution (n-k where n is the number of observations and k is the number of coefficients).

    Returns:
        pd.DataFrame: The dataframe with p-score and certainty columns.
    """
    # calculate p-scores and certainty
    predictions["p_score"] = predictions.apply(
        lambda row: p_score(row["pred"], row["se"], dof=dof, null_value=0),
        axis=1,
    )
    predictions["certainty"] = predictions["p_score"].apply(assign_certainty)
    return predictions


def extrapolate_predictions(df, year=2100):
    grouping_cols = ["scenario", "percentile", "core_grouping", "time_frame"]
    value_cols = [col for col in df.columns if col not in grouping_cols]

    new_rows = []

    for (scenario, percentile, core_grouping), group_df in df.groupby(
        ["scenario", "percentile", "core_grouping"]
    ):
        group_df = group_df[group_df["time_frame"] > 1995]

        if group_df.empty:
            continue

        interp_xs = group_df["time_frame"].values

        # prepare a dictionary for the new row (constant fields first)
        new_row = {
            "scenario": scenario,
            "percentile": percentile,
            "core_grouping": core_grouping,
            "time_frame": year,
        }

        for value_col in value_cols:
            inter_ys = group_df[value_col].values

            # need at least 2 points to interpolate/extrapolate
            if len(interp_xs) < 2:
                continue

            spline = make_interp_spline(
                interp_xs, inter_ys, k=min(2, len(interp_xs) - 1)
            )
            value_at_year = float(spline(year))  # returns as array

            new_row[value_col] = value_at_year

        new_rows.append(new_row)

    # add the new rows to the original dataframe
    df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

    return df


def fit_curve(
    df: pd.DataFrame, variable: str, effect_type: str, order: int
) -> sm.regression.linear_model.RegressionResultsWrapper:
    """
    Fit a polynomial curve to the data.

    Args:
        df (pd.DataFrame): The dataframe containing the data.
        variable (str): The independent variable.
        effect_type (str): The dependent variable.
        order (int): The order of the polynomial to fit.

    Returns:
        model: The fitted regression model.
    """
    df = df[df[variable].notna() & df[effect_type].notna()]

    # create polynomial features and fit model
    X = np.vander(df[variable], N=order + 1, increasing=True)
    model = sm.OLS(df[effect_type], X).fit()
    return model


def get_moderator_names(model: ro.vectors.ListVector) -> list[str]:
    """Get the names of the moderators from the model."""
    x = ro.r.assign("x", model)  # noqa
    beta_rownames = []
    for rn in ro.r("rownames(x$beta)"):
        beta_rownames.append(str(rn))
    return beta_rownames


def get_moderator_index(
    model: ro.vectors.ListVector, moderator_names: list[str]
) -> int:
    """Get the index of the moderator variable in the predictors list, accounting for intercept."""
    return (
        get_moderator_names(model).index(moderator_names)
        if isinstance(moderator_names, str)
        else [get_moderator_names(model).index(name) for name in moderator_names]
    )


def get_value_from_surface(model_surface, meshgrids, dt, dp):
    """Find the indices of the closest values to dt and dp in the sorted arrays"""
    i = np.abs(meshgrids[0][:, 0] - dt).argmin()
    j = np.abs(meshgrids[1][0, :] - dp).argmin()
    return model_surface[i, j]


def populate_anomaly_df_with_surface_values(
    anomaly_df, model_surface, meshgrids
) -> pd.DataFrame:
    for k in model_surface.keys():
        anomaly_df.loc[:, k] = anomaly_df.apply(
            lambda row: get_value_from_surface(
                model_surface[k],
                meshgrids,
                row["anomaly_value_sst"],
                row["anomaly_value_ph"],
            ),
            axis=1,
        )
    return anomaly_df
