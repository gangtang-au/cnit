# CNit

**CNit** (Carbon-Nitrogen Interactions in Terrestrial ecosystems) is a process-based terrestrial biogeochemistry model that simulates coupled carbon and nitrogen cycles in terrestrial ecosystems.  It is the land component of MAGICC (Model for the Assessment of Greenhouse Gas Induced Climate Change).

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
from cnit.physics.carbon_nitrogen_cycle import CarbonNitrogenCycleModel
from cnit.physics.carbon_nitrogen_cycle import CarbonNitrogenCycleModelConfig
from cnit.utils.units import Q

# Create model
config = CarbonNitrogenCycleModelConfig()
model = CarbonNitrogenCycleModel.from_config(config)

# Set up forcing
time_axis = Q(np.arange(1850, 2101), "yr")
dT_s = Q(np.linspace(0, 2, len(time_axis)), "K")
CO2_s = Q(np.linspace(280, 560, len(time_axis)), "ppm")
# ... other forcing inputs

# Run simulation
results = model.run(time_axis=time_axis, dT_s=dT_s, CO2_s=CO2_s, ...)

# Access results
print(f"Final plant carbon: {results['CplsP'][-1]}")
```

## Documentation
The full documentation for CNit is available at 
[https://cnit.readthedocs.io](https://cnit.readthedocs.io).

## License
BSD 3-Clause License. See LICENSE for details.

Copyright (c) 2025, Gang Tang and contributors.

## Authors and Contributors

**Author, Developer, and Maintainer:**  
Gang Tang  
- gang.tang.au@gmail.com  
- gang.tang@student.unimelb.edu.au  
- tgang@bgc-jena.mpg.de

**Contributors:**
- Zebedee Nicholls (zebedee.nicholls@climate-energy-college.org)
- Alexander Norton (alex.norton@csiro.au)
- Sönke Zaehle (szaehle@bgc-jena.mpg.de)
- Malte Meinshausen (malte.meinshausen@unimelb.edu.au)

## Citation
If you use CNit in your research, please cite:

Tang, G., Nicholls, Z., Norton, A., Zaehle, S., and Meinshausen, M.: Synthesizing
global carbon–nitrogen coupling effects – the MAGICC coupled carbon–nitrogen cycle
model v1.0, Geosci. Model Dev., 18, 2193–2230,
https://doi.org/10.5194/gmd-18-2193-2025, 2025.

## Related Projects

MAGICC - Model for the Assessment of Greenhouse Gas Induced Climate Change,
- https://magicc.org/
- https://gitlab.com/magicc/magicc
