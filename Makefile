.PHONY: maint test lint docs upload clean

maint:
	uv lock --upgrade
	uv run pre-commit autoupdate
	uv run pre-commit run --all-files

test:
	uv run pytest

lint:
	uv run ruff format --check
	uv run ruff check
	uv run mypy

docs:
	uv run --group docs sphinx-build -b html docs/source docs/build/html

upload:
	make clean
	uv build
	uv publish

clean:
	rm -rf aaindex1 aaindex2 aaindex3
	rm -rf build dist *.egg-info tests/reports docs/build .pytest_cache .mypy_cache .ruff_cache .coverage
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
