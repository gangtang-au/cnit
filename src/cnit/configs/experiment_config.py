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
    the config checks that ``time_axis`` has at least two finite, strictly increasing,
    evenly spaced coordinates and that every forcing series has exactly the same shape
    as the axis. It does not check the forcing-series values for finiteness or verify
    that unused net- and gross-land-use series are zero. The model's land-use method
    settings determine which series are used.

    """

    name: str
    """Label used to identify this scenario in the combined results dataset."""

    switch_N: pint.Quantity = field(validator=check_units("1"))
    """Dimensionless switch for nitrogen limitation and the nitrogen cycle: one
    enables coupling and zero disables it [dimensionless]."""

    time_axis: pint.Quantity = field(validator=check_units("yr"))
    """Evenly spaced time coordinates in years. The supplied calibrations use annual
    spacing, and evenly spaced subannual axes are also accepted [yr]."""

    dT_s: pint.Quantity = field(validator=check_units("K"))
    """Temperature anomaly relative to the model reference temperature [K]."""

    CO2_s: pint.Quantity = field(validator=check_units("ppm"))
    """Atmospheric CO₂ concentration [ppm]."""

    CemsLUnet_s: pint.Quantity = field(validator=check_units("GtC/yr"))
    """Net land-use carbon emission rate; used when the carbon land-use method
    selects net emissions [GtC/yr]."""

    CemsLUgrs_s: pint.Quantity = field(validator=check_units("GtC/yr"))
    """Gross land-use carbon emission rate; used when the carbon land-use method
    selects gross emissions [GtC/yr]."""

    NflxAD_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    """Nitrogen input rate from atmospheric deposition [GtN/yr]."""

    NflxFT_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    """Nitrogen input rate from fertilizer [GtN/yr]."""

    NemsLUnet_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    """Net land-use organic nitrogen emission rate; used when the nitrogen land-use
    method selects net emissions [GtN/yr]."""

    NemsLUgrs_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    """Gross land-use organic nitrogen emission rate; used when the nitrogen land-use
    method selects gross emissions [GtN/yr]."""

    NemsLUmin_s: pint.Quantity = field(validator=check_units("GtN/yr"))
    """Direct mineral nitrogen loss associated with land use [GtN/yr]."""

    metadata: Optional[Dict[str, object]] = field(default=None)
    """Optional user-provided metadata. It is stored on the config but is not passed
    as a forcing input to the model."""

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
