from typing import Dict, Optional

from attrs import define, field
import numpy as np
import pint

from ..utils.units import check_units


@define
class CarbonNitrogenCycleExperimentConfig:
    """Time axis and forcing series for one model experiment.

    Use this class with
    :meth:`~cnit.physics.carbon_nitrogen_cycle.CarbonNitrogenCycleModel.run_experiments`.
    The model configuration stores parameters shared by the model; each experiment
    config stores the time axis and forcing data for one scenario, along with its name
    and nitrogen switch.

    Field validators check that inputs have compatible physical units. At construction,
    the config also checks that ``time_axis`` has at least two finite, strictly
    increasing, evenly spaced coordinates and that every forcing series has exactly
    the same shape as the axis. It does not check forcing values for finiteness or
    verify that unused net- and gross-land-use series are zero. The model's land-use
    method settings determine which series are used.

    Attributes
    ----------
    name
        Label used to identify this scenario in the combined results dataset.
    switch_N
        Dimensionless switch for nitrogen limitation and the nitrogen cycle: one
        enables coupling and zero disables it.
    time_axis
        Time coordinates in years. Samples are evenly spaced; each forcing value
        represents the interval beginning at that coordinate. Annual spacing is used
        by the supplied calibrations, and evenly spaced subannual axes are accepted.
    dT_s
        Temperature anomaly relative to the model reference temperature [K].
    CO2_s
        Atmospheric CO₂ concentration [ppm].
    CemsLUnet_s, CemsLUgrs_s
        Net and gross land-use carbon emission rates [GtC/yr]. The model's carbon
        land-use method selects which series is read.
    NflxAD_s, NflxFT_s
        Nitrogen input rates from atmospheric deposition and fertilizer [GtN/yr].
    NemsLUnet_s, NemsLUgrs_s
        Net and gross land-use organic nitrogen emission rates [GtN/yr]. The model's
        nitrogen land-use method selects which series is read.
    NemsLUmin_s
        Direct mineral nitrogen loss associated with land use [GtN/yr].
    metadata
        Optional user-provided metadata. It is stored on the config but is not passed
        as a forcing input to the model.
    """

    name: str
    switch_N: pint.Quantity = field(validator=check_units("1"))
    time_axis: pint.Quantity = field(validator=check_units("yr"))
    dT_s: pint.Quantity = field(validator=check_units("K"))
    CO2_s: pint.Quantity = field(validator=check_units("ppm"))
    CemsLUnet_s: pint.Quantity = field(validator=check_units("GtC/yr"))
    CemsLUgrs_s: pint.Quantity = field(validator=check_units("GtC/yr"))
    NflxAD_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    NflxFT_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    NemsLUnet_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    NemsLUgrs_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    NemsLUmin_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    metadata: Optional[Dict[str, object]] = field(default=None)

    def __attrs_post_init__(self) -> None:
        """Check time-axis regularity and forcing-series lengths.

        Field validators check unit compatibility before this method runs.
        """
        time_values = np.asarray(self.time_axis.to("yr").m, dtype=float)
        if time_values.ndim != 1 or time_values.size < 2:
            raise ValueError(
                "time_axis must be one-dimensional with at least two points."
            )
        if not np.all(np.isfinite(time_values)):
            raise ValueError("time_axis must contain only finite values.")

        time_steps = np.diff(time_values)
        if np.any(time_steps <= 0):
            raise ValueError("time_axis must be strictly increasing.")
        if not np.allclose(time_steps, time_steps[0], rtol=1e-10, atol=1e-12):
            raise ValueError("time_axis must be evenly spaced.")

        series = {
            "dT_s": self.dT_s,
            "CO2_s": self.CO2_s,
            "CemsLUnet_s": self.CemsLUnet_s,
            "CemsLUgrs_s": self.CemsLUgrs_s,
            "NflxAD_s": self.NflxAD_s,
            "NflxFT_s": self.NflxFT_s,
            "NemsLUnet_s": self.NemsLUnet_s,
            "NemsLUgrs_s": self.NemsLUgrs_s,
            "NemsLUmin_s": self.NemsLUmin_s,
        }
        for name, values in series.items():
            if values.shape != self.time_axis.shape:
                raise ValueError(
                    f"{name} has shape {values.shape}, but time_axis has shape "
                    f"{self.time_axis.shape}."
                )

    def to_input_dict(self):
        """Return the model inputs stored in this experiment config.

        The returned dictionary can be expanded into ``CarbonNitrogenCycleModel.run``.
        The experiment ``name`` and ``metadata`` are not included.
        """
        return dict(
            time_axis=self.time_axis,
            dT_s=self.dT_s,
            CO2_s=self.CO2_s,
            CemsLUnet_s=self.CemsLUnet_s,
            CemsLUgrs_s=self.CemsLUgrs_s,
            NflxAD_s=self.NflxAD_s,
            NflxFT_s=self.NflxFT_s,
            NemsLUnet_s=self.NemsLUnet_s,
            NemsLUgrs_s=self.NemsLUgrs_s,
            NemsLUmin_s=self.NemsLUmin_s,
        )
