# Python makefile https://krzysztofzuraw.com/blog/2016/makefiles-in-python-projects.html
# Delcare all non-file targets as phony
.PHONY: clean clean-build clean-data data isort lint test
TEST_PATH=./

help:
	@echo "	clean-pyc"
	@echo "		Remove python artifacts."
	@echo "	clean-build"
	@echo "		Remove build artifacts."
	@echo "	clean-data"
	@echo "		Remove test data."
	@echo "	data"
	@echo "		Create test data."
	@echo "	atmos"
	@echo "		Prepare atmosphere model data."
	@echo "	isort"
	@echo "		Sort import statements."
	@echo "	lint"
	@echo "		Check style with flake8."
	@echo "	test"
	@echo "		Run py.test"
	@echo "	init"
	@echo "		Initialise by installing requirements"
	@echo "	init-dev"
	@echo "		Initialise by installing normal and dev requirements"
	@echo "	cov"
	@echo "		Produce coverage report"
	@echo "	mypy"
	@echo "		Run type checking with mypy"

clean-pyc:
	find . -name '*.pyc' -exec rm --force {} +
	find . -name '*.pyo' -exec rm --force {} +
	find . -name '*~' -exec rm --force  {} +

clean-data:
	rm --force --recursive data/test_data/PHOENIX-ACES_spectra
	rm --force --recursive data/test_data/results
	rm --force --recursive data/test_data/resampled

data:
	python eniric_scripts/make_test_data.py

atmos:
	split_atmmodel.py -b ALL
	bary_shift_atmmodel.py -b ALL

clean-build:
	rm --force --recursive build/
	rm --force --recursive dist/
	rm --force --recursive *.egg-info

isort:
	sh -c "isort  --recursive . "

lint:
	flake8 --exclude=.tox

test: clean-pyc
	py.test --verbose --color=yes $(TEST_PATH)

init:
	pip install -r requirements.txt

init-dev:
	pip install -r requirements.txt
	pip install -r requirements_dev.txt

cov: $(module)/*
	py.test --cov=$(module)
	coverage html

mypy:
	# py.test --mypy
	mypy --ignore-missing-imports .
