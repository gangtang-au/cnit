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

.. code-block:: bash

   pip install cnit

For Developers
~~~~~~~~~~~~~~

To set up a development environment for this project, follow these steps:

1. **Create a new conda environment** (replace ``cnit`` with your preferred name if desired):

   .. code-block:: bash

      conda create --name cnit python=3.9

2. **Activate your environment**:

   .. code-block:: bash

      conda activate cnit

3. **Clone the repository** (if you haven't already):

   .. code-block:: bash

      git clone https://github.com/gangtang-au/cnit.git
      cd cnit

4. **Install all dependencies and set up the package**:

   .. code-block:: bash

      make env

   This will install all required dependencies (including notebook, docs, and testing tools)
   and install the package in editable mode.

.. note::

   The Makefile will check that you are not in the ``base`` environment and will refuse
   to run if you are.

**TL;DR for developers:**

.. code-block:: bash

   conda create --name cnit python=3.9
   conda activate cnit
   git clone https://github.com/gangtang-au/cnit.git
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

The ``docs`` extra installs the Python documentation packages. Sphinx also needs
the external Pandoc executable to convert notebooks. Install Pandoc and ensure
``pandoc`` is on your ``PATH`` before building the documentation. For example:

.. code-block:: bash

   # macOS with Homebrew
   brew install pandoc

   # Ubuntu/Debian
   sudo apt-get install pandoc

   # Windows with winget
   winget install --source winget --exact --id JohnMacFarlane.Pandoc

See the `Pandoc installation guide <https://pandoc.org/installing.html>`_ for
other installation options. From a venv, build the documentation with:

.. code-block:: bash

   python scripts/build_docs.py

The repository-level ``make docs`` target requires an active Conda environment.

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

The runnable basic usage example is provided as a Jupyter notebook. It covers
model setup, illustrative forcing, running a simulation, inspecting results, and
plotting the carbon and nitrogen fluxes and pools.

.. toctree::
   :maxdepth: 1

   CNit_quickstart_basic_usage

Next Steps
----------

- Read :doc:`overview` to understand the model structure
- Browse :doc:`api/cnit.physics` for detailed API documentation
- Check :doc:`example` for CMIP6 ESM emulation and TRENDY DGVM emulation examples.
