from contextlib import contextmanager

# project_root = Path(__file__).parent.parent
# sys.path.insert(0, str(project_root))


@contextmanager
def RContextManager():
    """Safely initialise R context"""

    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter

        # Ensure pandas2ri is activated
        pandas2ri.activate()

        # Use the combined converter context
        with localconverter(ro.default_converter + pandas2ri.converter):
            yield {"ro": ro, "pandas2ri": pandas2ri, "localconverter": localconverter}
    except Exception as e:
        print(f"[R_CONTEXT] Error in R context: {e}")
        raise
