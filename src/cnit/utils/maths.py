"""
Mathematical functions and operations
"""
import numpy as np
import pint
from typing import Any, Callable
import scipy.interpolate

from .units import UNIT_REGISTRY


def calculate_cumulative_sum(
    in_arr: pint.Quantity,
    time_axis: pint.Quantity,
) -> pint.Quantity:
    """
    Calculate cumulative sum

    Parameters
    ----------
    in_arr
        Array of which to calculate the cumulative sum

    time_axis
        Points in time at which the points apply

    Returns
    -------
        Cumulative sum

    Notes
    -----
    The current implementation assumes that each value in ``in_arr`` applies
    for the entire timestep and is constant for the entire timestep. This is
    basically equivalent to a left-hand Riemann integral.

    We also assume that the last step is as long as the second to last step.
    Doing this better requires taking time bounds as input, not the time axis.

    TODO: This should be unified with the logic that is used in scmdata.
    """
    time_diff = time_axis[1:] - time_axis[:-1]
    time_diff = np.hstack([time_diff, time_diff[-1]])

    in_arr_total = in_arr * time_diff

    in_arr_cumulative_sum = np.cumsum(in_arr_total)

    return in_arr_cumulative_sum


def _get_interp(
    v: pint.Quantity,
    ta: pint.Quantity,
    **kwargs: Any,
) -> Callable[[pint.Quantity], pint.Quantity]:
    """
    Get interpolated function for pint quantities

    Parameters
    ----------
    v
        Array to interpolate

    ta
        Time axis

    **kwargs
        Passed to :func:`scipy.interpolate.interp1d`

    Returns
    -------
        Callable which gives interpolated values for
        arbitrary values of time
    """
    return UNIT_REGISTRY.wraps(v.units, ta.units)(
        scipy.interpolate.interp1d(
            ta.m, v.m, bounds_error=False, fill_value="extrapolate", **kwargs
        )
    )


class SolveError(ValueError):
    """
    Exception raised when a model can't be solved

    This could occur because the integration has a runaway effect in it
    """
