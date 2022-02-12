maint:
	pip install -r requirements/dev.txt
	pre-commit autoupdate && pre-commit run --all-files
	pip-compile -U setup.py
	pip-compile -U requirements/ci.in
	pip-compile -U requirements/dev.in

upload:
	make clean
	python setup.py sdist bdist_wheel && twine upload -s dist/*

clean:
	python setup.py clean --all
	pyclean .
	rm -rf aaindex1 aaindex2 aaindex3 propy/aaindex1 propy/aaindex2 propy/aaindex3
	rm -rf *.pyc __pycache__ build dist propy3.egg-info propy/__pycache__ tests/__pycache__ tests/reports docs/build .pytest_cache .tox .coverage

mypy:
	mypy . --ignore-missing-imports --python-version 3.8
