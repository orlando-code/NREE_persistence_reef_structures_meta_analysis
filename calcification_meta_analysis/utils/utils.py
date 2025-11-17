import logging
import string
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


def convert_png_to_jpg(directory: str | Path):
    """Convert all PNG files in a directory to JPG format.

    Args:
        directory (str | Path): The directory to convert PNG files to JPG format.

    Returns:
        generator: A generator of the converted JPG files.
    """
    png_files = Path(directory).glob("*.png")
    for png_file in png_files:
        with Image.open(png_file) as img:
            rgb_img = img.convert("RGB")
            jpg_file = png_file.with_suffix(".jpg")
            rgb_img.save(jpg_file, "JPEG")
            logging.info(f"Converted {png_file} to {jpg_file}")
            yield jpg_file


def uniquify_repeated_values(vals: list, uniquify_str: str = "LOC") -> list:
    """
    Append a unique suffix to repeated values in a list.

    Parameters:
        vals (list): List of values.

    Returns:
        list: List of values with unique suffixes.
    """
    from collections import Counter, defaultdict

    # precompute counts for all values to avoid repeated .count() calls
    value_counts = Counter(vals)
    counts = defaultdict(int)
    result = []
    for val in vals:
        count = counts[val]
        if value_counts[val] > 1:
            suffix = f"-{uniquify_str}-{string.ascii_uppercase[count]}"
        else:
            suffix = ""
        result.append(f"{val}{suffix}")
        counts[val] += 1
    return result


def safe_to_numeric(col):
    """Convert column to numeric if possible, otherwise return as is.

    Args:
        col (pd.Series): The column to convert to numeric.

    Returns:
        pd.Series: The converted column.
    """
    try:
        return pd.to_numeric(col)
    except (ValueError, TypeError):
        return col  # return original column if conversion fails


def aggregate_df(df: pd.DataFrame, method: str = "mean") -> pd.DataFrame:
    """Aggregate DataFrame by specified method (mean, median, etc.).

    Args:
        df (pd.DataFrame): The dataframe to aggregate.
        method (str): The method to use for aggregation (mean, median, etc.).

    Returns:
        pd.DataFrame: The aggregated dataframe.
    """
    aggregation_funcs = {
        col: method if pd.api.types.is_numeric_dtype(df[col]) else lambda x: x.iloc[0]
        for col in df.columns
    }  # define aggregation functions for each column
    return df.agg(aggregation_funcs)  # aggregate DataFrame


def calc_sd_from_se(se: float, n: int) -> float:
    """Calculate standard deviation from standard error and sample size

    Args:
        se (float): standard error
        n (int): number of samples

    Returns:
        float: standard deviation
    """
    return se * np.sqrt(n)


def round_down_to_nearest(x: float, step: float) -> float:
    """Round down to the nearest multiple of step.

    Args:
        x (float): The value to round down.
        step (float): The step size e.g. 0.5 for rounding to the nearest 0.5.
    Returns:
        float: The rounded down value.
    """
    return step * np.floor(x / step)


def get_unique_values(df: pd.DataFrame, column: str) -> list[str]:
    """Get unique values from a column.

    Args:
        df (pd.DataFrame): The dataframe to get unique values from.
        column (str): The column to get unique values from.

    Returns:
        list[str]: The unique values from the column. If no unique values are found, an empty list is returned.
    """
    if column in df.columns:
        unique_vals = df[column].dropna().unique().tolist()
        return sorted([str(v) for v in unique_vals])
    return []
