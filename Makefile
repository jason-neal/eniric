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
	@echo "	pdf"
	@echo "		Make pdf from paper.md"
	@echo "	publish_setup"
	@echo "		Update the tools for distribution."
	@echo "	build_dist"
	@echo "		Build the distribution of eniric."
	@echo "	pypi_test"
	@echo "		Publish eniric to Test PyPi."
	@echo "	pypi_publish"
	@echo "		Publish eniric to PyPi."


clean-pyc:
	find . -name '*.pyc' -exec rm --force {} +
	find . -name '*.pyo' -exec rm --force {} +
	find . -name '*~' -exec rm --force  {} +

clean-data:
	rm --force --recursive data/
	rm --force --recursive tests/data/phoenix-raw
	rm --force --recursive tests/data/btsettl-raw

atmos:
	split_atmmodel.py -b ALL
	barycenter_broaden_atmmodel.py -b ALL

clean-build:
	rm --force --recursive build/
	rm --force --recursive dist/
	rm --force --recursive *.egg-info

isort:
	sh -c "isort  --recursive . "

lint:
	flake8 --exclude=.tox

test: clean-pyc
	py.test --color=yes $(TEST_PATH)

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

pdf:
	pip install pandoc>=1.0.2
	(cd paper && pandoc --filter pandoc-citeproc --bibliography=paper.bib --variable classoption=onecolumn --variable papersize=a4paper -s paper.md -o paper.pdf)

publish_setup:
	pip install --upgrade twine setuptools wheel

build_dist:
	python setup.py sdist bdist_wheel

pypi_test:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

pypi_publish:
	twine upload dist/*
