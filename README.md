<table>
  <tr>
    <td valign="middle">
      <img src="https://raw.githubusercontent.com/gangtang-au/cnit/main/docs/_static/CNit_logo.png" alt="CNit logo" width="75">
    </td>
    <td valign="middle"><h1>CNit</h1></td>
  </tr>
</table>

[![CI](https://github.com/gangtang-au/cnit/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/gangtang-au/cnit/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/cnit)](https://pypi.org/project/cnit/)
[![Documentation](https://img.shields.io/readthedocs/cnit/latest)](https://cnit.readthedocs.io/en/latest/)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/cnit?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/cnit)

**CNit** (Carbon-Nitrogen Interactions in Terrestrial ecosystems) is a process-based
terrestrial biogeochemistry model of coupled carbon and nitrogen cycles. It is the land
component of **MAGICC** (Model for the Assessment of Greenhouse Gas Induced Climate
Change), designed to represent terrestrial carbon–nitrogen feedbacks efficiently in
integrated assessment and climate projection studies.

## Features

- Three carbon pools (plant, litter, and soil) and four nitrogen pools (plant, litter,
  soil, and mineral).
- Bidirectional carbon–nitrogen coupling: nitrogen availability limits carbon uptake,
  while carbon uptake drives nitrogen demand and can stimulate biological fixation.
- Effects of atmospheric CO₂, temperature, nitrogen availability, and land use on
  ecosystem fluxes.
- Separate deforestation, afforestation, post-deforestation regrowth, and afforestation
  decay fluxes for land-use and carbon dioxide removal scenarios.
- Annual timesteps in the default examples, with support for evenly spaced subannual
  timesteps. Allocation fractions apply within each model timestep; supplied calibration
  files are based on an annual timescale.

## Installation

Install CNit from PyPI:

```bash
pip install cnit
```

## Quick Start

This example runs CNit with simple annual input series. See the
[Quickstart Guide](https://cnit.readthedocs.io/en/latest/quickstart.html) for installation
options, a notebook version, and further explanation.

```python
import numpy as np
from cnit import CNitModel, CNitModelConfig, Q

# Create a model with default parameters.
model = CNitModel.from_config(CNitModelConfig())

# Set up annual time coordinates and input series.
time_axis = Q(np.arange(1850, 2101), "yr")
n_steps = len(time_axis)
dT_s = Q(np.linspace(0, 2, n_steps), "K")
CO2_s = Q(np.linspace(280, 560, n_steps), "ppm")

# Land-use emissions are zero for this simple example.
CemsLUnet_s = Q(np.zeros(n_steps), "GtC/yr")
CemsLUgrs_s = Q(np.zeros(n_steps), "GtC/yr")
NemsLUnet_s = Q(np.zeros(n_steps), "GtN/yr")
NemsLUgrs_s = Q(np.zeros(n_steps), "GtN/yr")
NemsLUmin_s = Q(np.zeros(n_steps), "GtN/yr")

# Nitrogen deposition and fertilizer inputs.
NflxAD_s = Q(np.linspace(0.05, 0.15, n_steps), "GtN/yr")
NflxFT_s = Q(np.linspace(0, 0.1, n_steps), "GtN/yr")

res = model.run(
    time_axis=time_axis,
    dT_s=dT_s,
    CO2_s=CO2_s,
    CemsLUnet_s=CemsLUnet_s,
    CemsLUgrs_s=CemsLUgrs_s,
    NflxAD_s=NflxAD_s,
    NflxFT_s=NflxFT_s,
    NemsLUnet_s=NemsLUnet_s,
    NemsLUgrs_s=NemsLUgrs_s,
    NemsLUmin_s=NemsLUmin_s,
)

print(f"Final plant carbon: {res['CplsP'].data[-1]} {res['CplsP'].units}")
print(f"Final plant nitrogen: {res['NplsP'].data[-1]} {res['NplsP'].units}")
print(f"Final NPP: {res['CflxNPP'].data[-1]} {res['CflxNPP'].units}")
```

## Documentation

The [CNit documentation](https://cnit.readthedocs.io/en/latest/) covers the model
overview, quickstart, detailed physics and configuration APIs, CMIP6 ESM and TRENDY DGVM
emulation workflows, and calibration data and results.

## References

If you use CNit in research, please cite the relevant publication(s):

### Model description and calibration

- **[CNit v1.0]** Tang, G., Nicholls, Z., Norton, A., Zaehle, S., and Meinshausen, M.:
  Synthesizing global carbon–nitrogen coupling effects – the MAGICC coupled
  carbon–nitrogen cycle model v1.0, *Geosci. Model Dev.*, 18, 2193–2230,
  [doi:10.5194/gmd-18-2193-2025](https://doi.org/10.5194/gmd-18-2193-2025), 2025.
- **[CNit v2.0]** Tang, G., Zaehle, S., Nicholls, Z., Norton, A., Ziehn, T., and
  Meinshausen, M.: Understanding the drivers of carbon–nitrogen cycle variability in
  CMIP6 ESMs with MAGICC CNit v2.0: model and calibration updates, *J. Adv. Model.
  Earth Syst.*, 18, e2025MS005270,
  [doi:10.1029/2025MS005270](https://doi.org/10.1029/2025MS005270), 2026.
- **[CNit v2.1]** Tang, G., Norton, A., Gasser, T., Nicholls, Z. R. J., Ziehn, T.,
  Abrahão, G. M., and Meinshausen, M. A.: Reassessing the terrestrial carbon sink
  reveals stronger CO₂ fertilization, elevated climate sensitivity, and intensifying
  nitrogen limitation, *ESS Open Archive*, preprint,
  [doi:10.22541/essoar.15006117/v2](https://doi.org/10.22541/essoar.15006117/v2), 2026.
- **[CMIP6 data]** Tang, G., Nicholls, Z., Jones, C., Gasser, T., Norton, A., Ziehn, T.,
  Romero-Prieto, A., and Meinshausen, M.: Investigating carbon and nitrogen conservation
  in reported CMIP6 Earth system model data, *Geosci. Model Dev.*, 18, 2111–2136,
  [doi:10.5194/gmd-18-2111-2025](https://doi.org/10.5194/gmd-18-2111-2025), 2025.

### Climate and carbon-cycle application

- **[MAGICC-CNit SSP projections]** Tang, G., Nicholls, Z., Norton, A., Gasser, T.,
  Ziehn, T., Zaehle, S., and Meinshausen, M.: Nitrogen limitation amplifies future
  warming by weakening terrestrial carbon cycle feedbacks and sink capacity,
  *Commun. Earth Environ.*, 7, 760,
  [doi:10.1038/s43247-026-03736-0](https://doi.org/10.1038/s43247-026-03736-0), 2026.

See [References](docs/project_info/references.rst) for the complete bibliography.

## Development

See [Contributing](docs/contributing.rst) for development setup and checks. The
[Changelog](docs/development/changelog.rst) lists release history.

## Authors and Contributors

**Author, developer, and maintainer:** Gang Tang

- [gang.tang.au@gmail.com](mailto:gang.tang.au@gmail.com) (primary)
- [gang.tang.1@unimelb.edu.au](mailto:gang.tang.1@unimelb.edu.au)
- [tgang@bgc-jena.mpg.de](mailto:tgang@bgc-jena.mpg.de)

**Contributors:**

- Zebedee Nicholls — [zebedee.nicholls@climate-energy-college.org](mailto:zebedee.nicholls@climate-energy-college.org)
- Alexander Norton — [alex.norton@csiro.au](mailto:alex.norton@csiro.au)
- Sönke Zaehle — [szaehle@bgc-jena.mpg.de](mailto:szaehle@bgc-jena.mpg.de)
- Malte Meinshausen — [malte.meinshausen@unimelb.edu.au](mailto:malte.meinshausen@unimelb.edu.au)

## Related Projects

CNit is the carbon cycle component of [MAGICC](https://magicc.org/). See the
[MAGICC repository](https://gitlab.com/magicc/magicc) for the climate model.

## License

CNit is distributed under the BSD 3-Clause License. See [LICENSE](LICENSE).
