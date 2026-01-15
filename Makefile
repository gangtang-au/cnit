.DEFAULT_GOAL := help

ifndef CONDA_PREFIX
$(error Conda not active, please install conda and then create an environment with `conda create --name cnit` and activate it with `conda activate cnit`)
else
ifeq ($(CONDA_DEFAULT_ENV),base)
$(error Do not install to conda base environment. Source a different conda environment e.g. `conda activate cnit` or `conda create --name cnit` and rerun make)
endif
VENV_DIR=$(CONDA_PREFIX)
endif

MAMBA_EXE := $(shell command -v mamba 2> /dev/null)
ifndef MAMBA_EXE
MAMBA_OR_CONDA=$(CONDA_EXE)
else
MAMBA_OR_CONDA=$(MAMBA_EXE)
endif

PYTHON=$(CONDA_PREFIX)/bin/python
CONDA_ENV_FILE=environment.yml
PYPROJECT_TOML=pyproject.toml

define PRINT_HELP_PYSCRIPT
import re, sys
for line in sys.stdin:
	match = re.match(r'^([\$$\(\)a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

help:  ## Print short description of each target
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

.PHONY: env
env:   ## Create/update conda environment and install package
	$(MAMBA_OR_CONDA) env update --name $(CONDA_DEFAULT_ENV) -f $(CONDA_ENV_FILE)
	$(VENV_DIR)/bin/pip install -e .[notebook,docs,tests]

.PHONY: docs
docs:  ## Build the docs
	$(MAKE) -C docs html

.PHONY: test
test:  ## Run the tests
	$(VENV_DIR)/bin/pytest tests -r a -vv

.PHONY: lint
lint:  ## Run flake8 and black for code style checks
	$(VENV_DIR)/bin/flake8 src/cnit
	$(VENV_DIR)/bin/black --check src/cnit tests

.PHONY: format
format:  ## Auto-format code with black
	$(VENV_DIR)/bin/black src/cnit tests

.PHONY: clean
clean:  ## Clean build, doc, and test artifacts
	rm -rf build/ dist/ .pytest_cache/ .mypy_cache/ .coverage htmlcov/ docs/_build/
	find . -type d -name '__pycache__' -exec rm -rf {} +