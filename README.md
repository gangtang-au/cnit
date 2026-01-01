<p align="left">
  <img src="docs/_static/CNit_logo.png" alt="CNit logo" width="100">
</p>

**CNit** (Carbon-Nitrogen Interactions in Terrestrial ecosystems) is a process-based 
terrestrial biogeochemistry model that simulates coupled carbon and nitrogen cycles 
in terrestrial ecosystems.  It is the land component of **MAGICC** (Model for the 
Assessment of Greenhouse Gas Induced Climate Change). 

## Features

* Explicit carbon-nitrogen coupling with nitrogen limitation feedbacks
* Four nitrogen pools (plant, litter, soil, mineral) and three carbon pools
* Environmental modifiers: CO₂ fertilization, temperature, land use
* Separation of deforestation and afforestation for CDR scenarios
* Annual timestep with sub-annual process representation

## Installation

**For users:**

```bash
pip install cnit
```

**For developers:**

To set up a development environment for this project, follow these steps:

1. **Create a new conda environment** (replace `cnit` with your preferred name if desired):

   ```bash
   conda create --name cnit python=3.9
   ```

2. **Activate your environment**:

   ```bash
   conda activate cnit
   ```

3. **Install all dependencies and set up the package**:

   ```bash
   make env
   ```

This will install all required dependencies (including notebook, docs, and testing tools) and install the package in editable mode.

> **Note:**  
> The Makefile will check that you are not in the `base` environment and will refuse to run if you are.

---

**TL;DR:**
```bash
conda create --name cnit python=3.9
conda activate cnit
make env
```

## Quick Start

```python
import numpy as np
from cnit import CNitModel, CNitModelConfig, CNitExpConfig
from cnit import Q # Local pint.quantity class for units handling

# Create model with default parameters
config = CNitModelConfig()
model = CNitModel.from_config(config)

# Set up time axis and forcing
time_axis = Q(np.arange(1850, 2101), "yr")
n_steps = len(time_axis)

# Simple forcing (constant values for illustration)
dT_s = Q(np.linspace(0, 2, n_steps), "K")  # 2°C warming
CO2_s = Q(np.linspace(280, 560, n_steps), "ppm")  # CO2 doubling

# Land use emissions (zero for simplicity)
CemsLUnet_s = Q(np.zeros(n_steps), "GtC/yr")
CemsLUgrs_s = Q(np.zeros(n_steps), "GtC/yr")

# Nitrogen inputs
NflxAD_s = Q(np.linspace(0.05, 0.15, n_steps), "GtN/yr")  # Increasing deposition
NflxFT_s = Q(np.linspace(0, 0.1, n_steps), "GtN/yr")  # Increasing fertilizer

# Nitrogen land use emissions
NemsLUnet_s = Q(np.zeros(n_steps), "GtN/yr")
NemsLUgrs_s = Q(np.zeros(n_steps), "GtN/yr")
NemsLUmin_s = Q(np.zeros(n_steps), "GtN/yr")

# Run the model
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

# Access results
cveg, cveg_unit = res["CplsP"].data, res["CplsP"].units
nveg, nveg_unit = res["NplsP"].data, res["NplsP"].units
npp, npp_unit = res["CflxNPP"].data, res["CflxNPP"].units

print(f"Final plant carbon: {cveg[-1]} {cveg_unit}")
print(f"Final plant nitrogen: {nveg[-1]} {nveg_unit}")
print(f"NPP change: {npp[0]} → {npp[-1]} {npp_unit}")
```

## Documentation
The full documentation for CNit is available at 
[https://cnit.readthedocs.io](https://cnit.readthedocs.io).

## License
BSD 3-Clause License. See LICENSE for details.

Copyright (c) 2026, Gang Tang and contributors.

## Authors and Contributors

**Author, Developer, and Maintainer:**  
Gang Tang  
- gang.tang.au@gmail.com (primary)
- gang.tang@student.unimelb.edu.au  
- tgang@bgc-jena.mpg.de

**Contributors:**
- Zebedee Nicholls (zebedee.nicholls@climate-energy-college.org)
- Alexander Norton (alex.norton@csiro.au)
- Sönke Zaehle (szaehle@bgc-jena.mpg.de)
- Malte Meinshausen (malte.meinshausen@unimelb.edu.au)

## Citation
If you use CNit in your research, please cite:

- Tang, G., Nicholls, Z., Norton, A., Zaehle, S., and Meinshausen, M.: Synthesizing
global carbon–nitrogen coupling effects – the MAGICC coupled carbon–nitrogen cycle
model v1.0, Geosci. Model Dev., 18, 2193–2230,
https://doi.org/10.5194/gmd-18-2193-2025, 2025.

- Tang, G., Zaehle, S., Nicholls, Z., Norton, A., Ziehn, T., & Meinshausen, 
M. Understanding the Drivers of Carbon-Nitrogen Cycle Variability in CMIP6 ESMs 
with MAGICC CNit v2.0: Model and Calibration Updates. ESS Open Archive. 
June 16, 2025. https://doi.org/10.22541/essoar.175008280.09297369/v1 
[accepted by Journal of Advances in Modeling Earth Systems (JAMES)]

## Related Projects

MAGICC - Model for the Assessment of Greenhouse Gas Induced Climate Change,
- https://magicc.org/
- https://gitlab.com/magicc/magicc
