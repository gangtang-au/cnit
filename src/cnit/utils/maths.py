"""
Mathematical functions and operations
"""

from typing import Any, Callable, List, Union

import numpy as np
import pint
import scipy.interpolate


class PackedLinearInterpolator:
    """Fast vectorized linear interpolator for multiple synchronized time series.

    Supports 1D arrays (single series) and 2D arrays (shape [n_series, n_times]).
    Supports both evenly spaced grids (O(1) index calculation) and non-uniform grids.
    """

    def __init__(self, time_values: np.ndarray, data: np.ndarray):
        self.t_eval = np.asarray(time_values, dtype=float)
        self.t0 = self.t_eval[0]
        self.n_times = len(self.t_eval)
        if self.n_times < 2:
            raise ValueError(
                "PackedLinearInterpolator requires at least two time points."
            )
        self.dt = (self.t_eval[-1] - self.t_eval[0]) / (self.n_times - 1)
        self.is_uniform = np.allclose(
            np.diff(self.t_eval), self.dt, rtol=1e-5, atol=1e-8
        )
        self.data = np.asarray(data, dtype=float)

    def __call__(self, t: float) -> Union[float, np.ndarray]:
        if self.is_uniform:
            i = int((t - self.t0) / self.dt)
        else:
            i = np.searchsorted(self.t_eval, t, side="right") - 1

        if i < 0:
            i = 0
        elif i >= self.n_times - 1:
            i = self.n_times - 2

        dt_i = self.t_eval[i + 1] - self.t_eval[i]
        alpha = (t - self.t_eval[i]) / dt_i
        return (1.0 - alpha) * self.data[..., i] + alpha * self.data[..., i + 1]


def _make_packed_numeric_interp(
    series_list: List[Any],
    time_values: np.ndarray,
    units: List[str],
) -> PackedLinearInterpolator:
    """Pack multiple quantity/numeric series into a single PackedLinearInterpolator."""
    data = []
    for s, u in zip(series_list, units):
        if isinstance(s, pint.Quantity):
            data.append(s.to(u).m)
        else:
            data.append(np.asarray(s, dtype=float))
    return PackedLinearInterpolator(time_values, np.vstack(data))


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
