Contributing
============

CNit is developed and maintained by Gang Tang as part of the MAGICC project.

How to Contribute
-----------------

We welcome bug reports, feature suggestions, and pull requests. For substantial
changes, open an issue first so the proposed work can be discussed before you
begin.

To report a bug or suggest a feature:

1. **Report bugs:** Open an issue on `GitHub <https://github.com/gangtang-au/cnit/issues>`_
2. **Suggest features:** Describe your idea in a GitHub issue
3. **Ask questions:** Email gang.tang.au@gmail.com or open a discussion

For a pull request, keep the change focused, add or update tests for changed
behavior, and update the documentation when user-facing behavior or examples
change. Include a short description of the change and the checks you ran.

Development Setup
-----------------

The commands below use Conda and the repository's ``environment.yml``. From the
repository root, run:

.. code-block:: bash

   git clone https://github.com/gangtang-au/cnit.git
   cd cnit

   # Create and activate a Python 3.9–3.12 environment
   conda create --name cnit python=3.12
   conda activate cnit

   # Install CNit and all development dependencies
   make env

``make env`` requires an active, non-``base`` Conda environment. It installs
CNit in editable mode with the notebook, documentation, test, lint, and format extras.
Python package dependencies and extras are defined in ``pyproject.toml``;
``environment.yml`` adds Conda-managed tools such as Pandoc.

Development checks
------------------

Run the checks relevant to your change before opening a pull request:

.. code-block:: bash

   make test
   make lint

   # Check formatting without changing files
   black --check src/cnit tests

   # Apply Black formatting when needed
   make format

   # Build the Sphinx documentation (requires Pandoc for notebook conversion)
   make docs

The test target runs the full test suite. The README Quick Start example is
executed as part of the regular tests to catch example drift. The documentation
build renders the notebooks with execution disabled, so notebook outputs should
be saved in the notebook files.

The model-results regression test compares the default model configuration with
the local ``main`` branch using the ``1pctCO2`` and ``hist_ssp585`` forcing
experiments. It also runs a calibrated ACCESS-ESM1-5 configuration and reports
maximum absolute and relative differences for selected carbon and nitrogen
diagnostics. To show those diagnostic values locally, run:

.. code-block:: bash

   CNIT_BASELINE_REF=origin/main python -m pytest tests/test_model_results_regression.py -q -s

Set ``CNIT_BASELINE_REF`` to another available Git ref when comparing against a
different baseline.

The repository-level ``notebooks/`` directory contains the notebook sources
used by the documentation build, which stages them in a temporary Sphinx
source directory.

Pull request checklist
----------------------

* Keep generated documentation output and local environment files out of commits.
* Run the tests, formatting check, lint, and relevant documentation checks.
* Describe compatibility or scientific-result changes clearly in the pull
  request.

CNit is released under the BSD 3-Clause License. See the
:doc:`License <../project_info/license>` page or the repository's LICENSE file
for details.
