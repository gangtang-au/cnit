"""
Mathematical functions and operations
"""

from typing import Any, Callable, Union

import numpy as np
import pint
import scipy.interpolate


def _make_numeric_interp(
    values: Any,
    time_values: np.ndarray,
    unit: str,
    *,
    extrapolate: bool,
) -> Callable[[Union[float, np.ndarray]], Union[float, np.ndarray]]:
    """Build a unit-free, linearly interpolated function from a quantity series."""
    if isinstance(values, pint.Quantity):
        values = values.to(unit).m
    else:
        values = np.asarray(values)

    return scipy.interpolate.interp1d(
        time_values,
        values,
        bounds_error=not extrapolate,
        fill_value="extrapolate" if extrapolate else np.nan,
    )


class SolveError(ValueError):
    """
    Exception raised when a model can't be solved

    This could occur because the integration has a runaway effect in it
    """
