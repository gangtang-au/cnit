from dataclasses import dataclass
from typing import Optional, Dict
import pint


@dataclass(slots=True)
class CarbonNitrogenCycleExperimentConfig:
    """
    Configuration for a Carbon-Nitrogen Cycle experiment.

    Attributes:
        name: Name of the experiment.
        time_axis: Time axis for the experiment [years].
        dT_s: Temperature anomaly time series [K].
        CO2_s: Atmospheric CO2 concentration time series [ppm].
        CemsLUnet_s: Net carbon emissions from land use [PgC/yr].
        CemsLUgrs_s: Gross carbon emissions from land use [PgC/yr].
        NflxAD_s: Nitrogen flux from atmospheric deposition [TgN/yr].
        NflxFT_s: Nitrogen flux from fertilizer [TgN/yr].
        NemsLUnet_s: Net nitrogen emissions from land use [TgN/yr].
        NemsLUgrs_s: Gross nitrogen emissions from land use [TgN/yr].
        NemsLUmin_s: Nitrogen emissions from mineralization [TgN/yr].
        metadata: Optional dictionary for additional metadata.
    """
    name: str
    switch_N: pint.Quantity
    time_axis: pint.Quantity
    dT_s: pint.Quantity
    CO2_s: pint.Quantity
    CemsLUnet_s: pint.Quantity
    CemsLUgrs_s: pint.Quantity
    NflxAD_s: pint.Quantity
    NflxFT_s: pint.Quantity
    NemsLUnet_s: pint.Quantity
    NemsLUgrs_s: pint.Quantity
    NemsLUmin_s: pint.Quantity
    metadata: Optional[Dict[str, object]] = None

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
