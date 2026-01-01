Quickstart Guide
================

This guide will get you up and running with CNit in minutes.

Installation
------------

Requirements
~~~~~~~~~~~~

- Python 3.9 or later
- pip (Python package installer)
- For developers: conda or mamba (recommended)

For Users
~~~~~~~~~

If you just want to use CNit:

..  code-block:: bash

   pip install cnit

For Developers
~~~~~~~~~~~~~~

To set up a development environment for this project, follow these steps:

1. **Create a new conda environment** (replace ``cnit`` with your preferred name if desired):

   .. code-block:: bash

      conda create --name cnit python=3.9

2. **Activate your environment**:

   .. code-block::  bash

      conda activate cnit

3. **Clone the repository** (if you haven't already):

   .. code-block:: bash

      git clone https://github.com/yourusername/cnit.git
      cd cnit

4. **Install all dependencies and set up the package**:

   .. code-block:: bash

      make env

   This will install all required dependencies (including notebook, docs, and testing tools)
   and install the package in editable mode.

..  note::

   The Makefile will check that you are not in the ``base`` environment and will refuse
   to run if you are.

**TL;DR for developers:**

..  code-block:: bash

   conda create --name cnit python=3.9
   conda activate cnit
   git clone https://github.com/yourusername/cnit.git
   cd cnit
   make env

Alternative: Using venv
~~~~~~~~~~~~~~~~~~~~~~~

If you prefer not to use conda:

.. code-block:: bash

   # Create a new virtual environment
   python -m venv cnit-env

   # Activate it
   # On macOS/Linux:
   source cnit-env/bin/activate
   # On Windows:
   cnit-env\Scripts\activate

   # Install in development mode
   pip install -e ".[notebook,docs,tests]"

Verify Installation
~~~~~~~~~~~~~~~~~~~

Check that CNit is installed correctly:

.. code-block:: bash

   pip show cnit

This will display the package information including version, location, and dependencies.

You can also verify by importing the package in Python:

.. code-block:: python

   import cnit
   print("CNit installed successfully!")

If the import succeeds without errors, CNit is ready to use!

Basic Usage
-----------

Import the main classes:

.. code-block::  python

   from cnit import CNModel, CNModelConfig, CNExpConfig

..  note::

   **Import Aliases:** For convenience, CNit provides short aliases:

   * ``CNModel`` = ``CarbonNitrogenCycleModel``
   * ``CNModelConfig`` = ``CarbonNitrogenCycleModelConfig``
   * ``CNExpConfig`` = ``CarbonNitrogenCycleExperimentConfig``

   Both names work identically.  Use whichever you prefer!

Your First Simulation
---------------------

Here's a minimal example to run CNit:

.. code-block::  python

   import numpy as np
   from cnit import CNitModel, CNitModelConfig
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
   cveg, cveg_unit = res["CplsP"].data, re["CplsP"].units
   nveg, nveg_unit = res["NplsP"].data, res["NplsP"].units
   npp, npp_unit = res["CflxNPP"].data, res["CflxNPP"].units

   print(f"Final plant carbon: {cveg[-1]} {cveg_unit}")
   print(f"Final plant nitrogen: {nveg[-1]} {nveg_unit}")
   print(f"NPP change: {npp[0]} → {npp[-1]} {npp_unit}")

Visualizing Results
-------------------

Plot the carbon pools over time:

.. code-block:: python

   import matplotlib.pyplot as plt

    # Set default plotting properties
    plt.rcParams.update(
        {
            "font.size": 6,  # Default font size for all text
            "axes.titlesize": 6,  # Axes title
            "axes.labelsize": 6,  # Axes labels
            "xtick.labelsize": 6,  # X-axis tick labels
            "ytick.labelsize": 6,  # Y-axis tick labels
            "legend.fontsize": 6,  # Legend text
            "figure.titlesize": 6,  # Figure title
            "lines.linewidth": 0.5,  # Default line width
            "figure.figsize": (9 / 2.54, 6 / 2.54),  # Default figure size in cm
            "figure.dpi": 150,  # Default DPI for figures
            "figure.constrained_layout.use": True,  # Use constrained layout by default
            "xtick.major.size": 1,  # X-axis major tick length
            "ytick.major.size": 1,  # Y-axis major tick length
        }
    )


    fig, (axs1, axs2) = plt.subplots(2, 4, figsize=(12 / 2.54, 6 / 2.54))
    time = time_axis.m
    # Carbon pools
    axs1[0].plot(time, res["CplsP"])
    axs1[1].plot(time, res["CplsL"])
    axs1[2].plot(time, res["CplsS"])
    axs1[3].axis("off")
    axs1[0].set_ylabel("Carbon Pool (GtC)")
    axs1[0].set_title("Plant")
    axs1[1].set_title("Litter")
    axs1[2].set_title("Soil")

    # Nitrogen pools
    axs2[0].plot(time, res["NplsP"])
    axs2[1].plot(time, res["NplsL"])
    axs2[2].plot(time, res["NplsS"])
    axs2[3].plot(time, res["NplsM"])
    axs2[0].set_ylabel("Nitrogen Pool (GtN)")
    axs2[0].set_title("Plant")
    axs2[1].set_title("Litter")
    axs2[2].set_title("Soil")
    axs2[3].set_title("Mineral")

    plt.show()

Understanding the Results
-------------------------

The results Dataset contains:

**Carbon Variables:**

- ``CplsP``, ``CplsL``, ``CplsS``: Plant, litter, and soil carbon pool sizes [GtC]
- ``CflxNPP``, ``CflxLPR``: Net primary production and litter production respiration [GtC/yr]
- ``CflxRH``: Heterotrophic respiration [GtC/yr]
- ``CflxNetPLS``: Net terrestrial carbon flux (net biome production) [GtC/yr]

**Nitrogen Variables:**

- ``NplsP``, ``NplsL``, ``NplsS``, ``NplsM``: Plant, litter, soil, and mineral nitrogen pool sizes [GtN]
- ``NflxPU``, ``NflxBNF``: Nitrogen plant uptake and biological nitrogen fixation [GtN/yr]
- ``NflxNetMIN``: Net mineralization [GtN/yr]
- ``NflxPUdef``: Nitrogen deficit for plant uptake [GtN/yr]

Next Steps
----------

- See :doc:`examples/index` for more detailed examples
- Read :doc:`overview` to understand the model structure
- Browse :doc:`api/cnit.physics` for detailed API documentation
- Check :doc:`example` for realistic forcing data