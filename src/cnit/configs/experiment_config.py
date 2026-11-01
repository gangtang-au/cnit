"""Configuration of the time axis and input data for CNit experiments.

The experiment input reference and attribute documentation are in
:class:`CarbonNitrogenCycleExperimentConfig`.
"""

from typing import Dict, Optional

from attrs import define, field
import numpy as np
import pint

from ..utils.units import check_units


@define
class CarbonNitrogenCycleExperimentConfig:
    """Stores the time axis and input series for one named CNit experiment.

    Use this config with
    :meth:`~cnit.physics.carbon_nitrogen_cycle.CarbonNitrogenCycleModel.run_experiments`.
    The model config stores parameters shared across runs; this config stores the name,
    time axis, input series, and nitrogen-cycle switch for one scenario. The public
    alias is ``CNitExpConfig``, available from the package root::

        from cnit import CNitExpConfig

    Field validators check unit compatibility. On construction, the config checks that
    ``time_axis`` has at least two finite, strictly increasing, evenly spaced
    coordinates and that every input series has exactly the same shape as the axis
    and contains only finite values. This includes unused input series: supply finite
    placeholders (for example, zeros) for those inputs. Unused net- and gross-land-use
    series need not be zero; the model's land-use method settings determine which
    series are used.

    Annual time axes are used in the supplied examples and calibration files. Evenly
    spaced subannual axes are also supported, and allocation fractions apply within
    each model timestep. Net and gross land-use series are both required, although the
    model uses only the series selected by the corresponding land-use method.

    The table summarizes the experiment inputs; the class attributes below provide
    their detailed documentation.

    .. list-table:: Experiment inputs
       :class: config-table
       :header-rows: 1
       :widths: 28 17 55

       * - Input
         - Unit
         - Description
       * - ``name``
         - —
         - Scenario label stored in the experiment results.
       * - ``switch_N``
         - dimensionless
         - Enables nitrogen-cycle coupling for this experiment when 1; disables it when 0.
       * - ``time_axis``
         - yr
         - Evenly spaced time coordinates.
       * - ``dT_s``
         - K
         - Temperature anomaly relative to the model reference temperature.
       * - ``CO2_s``
         - ppm
         - Atmospheric CO₂ concentration.
       * - ``CemsLUnet_s``
         - GtC/yr
         - Net carbon land-use emission rate; used when ``method_cLU`` is 0.
       * - ``CemsLUgrs_s``
         - GtC/yr
         - Gross carbon land-use emission rate; used when ``method_cLU`` is 1.
       * - ``NflxAD_s``
         - GtN/yr
         - Nitrogen input rate from atmospheric deposition.
       * - ``NflxFT_s``
         - GtN/yr
         - Nitrogen input rate from fertilizer.
       * - ``NemsLUnet_s``
         - GtN/yr
         - Net organic nitrogen land-use emission rate; used when ``method_nLU`` is 0.
       * - ``NemsLUgrs_s``
         - GtN/yr
         - Gross organic nitrogen land-use emission rate; used when ``method_nLU`` is 1.
       * - ``NemsLUmin_s``
         - GtN/yr
         - Direct mineral nitrogen loss associated with land use.
       * - ``metadata``
         - —
         - Optional user metadata stored on the config and omitted from model inputs.
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
    as an input to the model."""

    def __attrs_post_init__(self) -> None:
        """Check time-axis regularity, input-series lengths and finite values.

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
            if not np.all(np.isfinite(values.m)):
                raise ValueError(f"{name} must contain only finite values.")

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
