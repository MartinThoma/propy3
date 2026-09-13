## 2.0.1

* MAINT: Publish releases to PyPI from GitHub Actions via Trusted Publishing.

## 2.0.0

Breaking changes:

* Python 3.10 or newer is required (Python 3.6 - 3.9 are no longer supported).
* `AAIndex.init()` / `AAIndex.get()` without a path now read the aaindex files
  bundled with propy instead of trying to download them via FTP.

Other changes:

* BUG: propy no longer depends on `pkg_resources`, which is missing on Python 3.12+
  environments without setuptools and removed from recent setuptools releases.
  Package data is now loaded with `importlib.resources`.
* MAINT: Downloads from UniProt and AAindex now use HTTPS.
* MAINT: Support Python 3.10 - 3.14.
* MAINT: Move packaging metadata to `pyproject.toml`; use uv, ruff and mypy
  for development; replace Travis CI with GitHub Actions.

## 1.1.1

* BUG: Fix Grantham data (#22)
* DOC: Installation instructions for conda

## 1.1.0

* Feature: PyPro.GetALL now has parameters
* BUG: GetAPseudoAAC2 used PAAC and now uses APAAC

## 1.0.2

* BUG: The CTD._Polarity and CTD._NormalizedVDWV values were wrong (see [change](https://github.com/MartinThoma/propy3/commit/6788e96c4aed77dad52d7dafce447c522d15b012)). Kudos to [Caio Fontes](https://github.com/Caiofcas) for reporting the issue, [Qiaole He](https://github.com/KimHe) for providing a fix, and [Jonathan Chen](https://github.com/jowch) for confirming + adding a PR to fix it.

## 1.0.1

* BUG: Fix Mutability.json (K was -56, but needs to be 56)
