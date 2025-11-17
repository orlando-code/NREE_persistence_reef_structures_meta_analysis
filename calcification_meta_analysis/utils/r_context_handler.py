from contextlib import contextmanager

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr


def safe_import_r_package(package_name: str, install_if_missing: bool = True):
    """
    Safely import an R package with clean error handling.

    Args:
        package_name (str): Name of the R package to import
        install_if_missing (bool): Whether to install the package if it's missing

    Returns:
        The imported R package object

    Raises:
        ImportError: If package cannot be installed or imported
    """
    with RContextManager():
        # with RContextManager() as r_ctx:
        # ro = r_ctx["ro"]

        try:
            # try to import the package within the context
            return importr(package_name)
        except Exception as e:
            # if import fails, try installation if requested
            if install_if_missing:
                print(f"Installing R package '{package_name}'...")

                try:
                    # use R's utils package to install from CRAN
                    utils = importr("utils")
                    utils.chooseCRANmirror(ind=1)  # Choose first CRAN mirror
                    utils.install_packages(package_name)
                    print(f"Successfully installed '{package_name}'")

                    # try importing again after installation
                    return importr(package_name)
                except Exception as install_error:
                    raise ImportError(
                        f"Failed to install and import '{package_name}': {install_error}"
                    )
            else:
                raise ImportError(f"Failed to import '{package_name}': {e}")


@contextmanager
def RContextManager():
    """
    Safely initialize R context with robust error handling.
    This function handles different versions of rpy2 and provides fallback mechanisms.
    """
    try:
        # try the newer rpy2 API first (3.5+)
        try:
            from rpy2.robjects.conversion import localconverter

            # use the combined converter context
            with localconverter(ro.default_converter + pandas2ri.converter):
                yield {
                    "ro": ro,
                    "pandas2ri": pandas2ri,
                    "localconverter": localconverter,
                }

        except ImportError:
            # fallback for older rpy2 versions
            try:
                from rpy2.robjects.conversion import Converter

                # create converter instance
                converter_instance = Converter(
                    ro.default_converter + pandas2ri.converter
                )

                # use the combined converter context
                with converter_instance.context():
                    yield {"ro": ro, "pandas2ri": pandas2ri, "Converter": Converter}

            except Exception:
                # final fallback - just activate pandas2ri without context management
                yield {"ro": ro, "pandas2ri": pandas2ri, "localconverter": None}

    except Exception as e:
        print(f"[R_CONTEXT] Error in R context: {e}")
        print(
            f"[R_CONTEXT] rpy2 version: {ro.__version__ if hasattr(ro, '__version__') else 'unknown'}"
        )
        raise


def index_named_list(named_list, name: str):
    """Index an R NamedList object by name of variable. Necessary since OrdDict objects now replaced by NamedList objects."""
    names = named_list._NamedList__names
    return named_list[names.index(name)]
