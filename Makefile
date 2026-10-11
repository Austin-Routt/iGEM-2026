# UH iGEM 2026 — AKGut flux-balance analysis
.PHONY: help install figures verify clean

help:
	@echo "make install   - install pinned dependencies into the active environment"
	@echo "make figures   - regenerate all six wiki figures into figures/"
	@echo "make verify    - check the model loads and headline fluxes reproduce"
	@echo "make clean     - remove generated figures"

install:
	pip install -r requirements.txt

figures:
	python regenerate_model_figures.py

verify:
	python scripts/verify_iHM1533_numbers.py

clean:
	rm -f figures/fig*.png figures/fig*.svg
