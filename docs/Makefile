# Keep the traditional docs-local command routed through the canonical build script.
ROOT_DIR := $(abspath ..)
PYTHON ?= python

.PHONY: html
html:
	$(PYTHON) "$(ROOT_DIR)/scripts/build_docs.py"
