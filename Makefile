upload:
	make clean
	python3 setup.py sdist bdist_wheel && twine upload dist/*

clean:
	python setup.py clean --all
	pyclean .
	rm -rf aaindex1 aaindex2 aaindex3
	rm -rf *.pyc __pycache__ build dist propy3.egg-info propy/__pycache__ tests/__pycache__ tests/reports docs/build .pytest_cache .tox .coverage

mypy:
	mypy . --ignore-missing-imports --python-version 3.8
