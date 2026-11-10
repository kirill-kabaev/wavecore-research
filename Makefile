
PAPER_DIR := papers/01_maxwell_helmholtz

.PHONY: paper clean check

paper:
	cd $(PAPER_DIR) && \
	latexmk -lualatex \
	-interaction=nonstopmode \
	-halt-on-error main.tex

clean:
	cd $(PAPER_DIR) && latexmk -C

check:
	git status --short
