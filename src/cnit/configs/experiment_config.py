from typing import Dict, Optional

from attrs import define, field
import numpy as np
import pint

from ..utils.units import check_units


@define
class CarbonNitrogenCycleExperimentConfig:
    """
    Configuration for a Carbon-Nitrogen Cycle experiment.

    Field validators check that each input has compatible units. When the instance
    is created, ``__attrs_post_init__`` checks that ``time_axis`` is finite,
    strictly increasing, and evenly spaced, and that every forcing series has the
    same number of values as the axis. It does not inspect forcing values or require
    a particular array shape beyond matching sizes.

    Attributes:
        name: Name of the experiment.
        time_axis: Time axis for the experiment [years].
        dT_s: Temperature anomaly time series [K].
        CO2_s: Atmospheric CO2 concentration time series [ppm].
        CemsLUnet_s: Net carbon emissions from land use [GtC/yr].
        CemsLUgrs_s: Gross carbon emissions from land use [GtC/yr].
        NflxAD_s: Nitrogen flux from atmospheric deposition [GtN/yr].
        NflxFT_s: Nitrogen flux from fertilizer [GtN/yr].
        NemsLUnet_s: Net nitrogen emissions from land use [GtN/yr].
        NemsLUgrs_s: Gross nitrogen emissions from land use [GtN/yr].
        NemsLUmin_s: Nitrogen emissions from mineralization [GtN/yr].
        metadata: Optional dictionary for additional metadata.
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
