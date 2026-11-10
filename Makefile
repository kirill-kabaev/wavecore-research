
SHELL := /bin/bash

PAPER_DIR := papers/01_maxwell_helmholtz
PYTHON := .venv/bin/python

.DEFAULT_GOAL := help

.PHONY: help paper clean ci site serve

help:
	@echo "WaveCore Research"
	@echo "  make paper  - Build LaTeX PDF"
	@echo "  make ci     - Run full local CI"
	@echo "  make site   - Build MkDocs website"
	@echo "  make serve  - Start local website"
	@echo "  make clean  - Clean LaTeX files"

paper:
	cd $(PAPER_DIR) && latexmk -lualatex \
		-interaction=nonstopmode \
		-file-line-error \
		-halt-on-error main.tex

ci:
	PATH="$(CURDIR)/.venv/bin:$$PATH" bash scripts/local_ci.sh

site:
	$(PYTHON) -m mkdocs build --strict

serve:
	$(PYTHON) -m mkdocs serve

clean:
	cd $(PAPER_DIR) && latexmk -C
