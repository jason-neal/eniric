# Python makefile https://krzysztofzuraw.com/blog/2016/makefiles-in-python-projects.html
# Delcare all non-file targets as phony
.PHONY: clean clean-build clean-data data isort lint test
TEST_PATH=./

help:
	@echo "	clean"
	@echo "		Remove python artifacts."
	@echo "	clean-build"
	@echo "		Remove build artifacts."
	@echo "	clean-data"
	@echo "		Remove test data."
	@echo "	data"
	@echo "		Create test data."
	@echo "	isort"
	@echo "		Sort import statements."
	@echo "	lint"
	@echo "		Check style with flake8."
	@echo "	test"
	@echo "		Run py.test"

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
