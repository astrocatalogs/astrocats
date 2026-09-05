# Astrocats: Open Astronomy Catalogs #

[![CI](https://github.com/astrocatalogs/astrocats/actions/workflows/ci.yml/badge.svg)](https://github.com/astrocatalogs/astrocats/actions/workflows/ci.yml)
[![Python Version](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org)
[![arXiv](https://img.shields.io/badge/arXiv-1605.01054-green.svg?style=flat)](http://arxiv.org/abs/1605.01054)

The Astrocats package enables astronomers to construct their own curated catalogs of astronomical data with the intention of producing shareable catalogs of that data in human-readable formats. Astrocats is used by several existing open astronomy catalogs, including:

* [The Open Supernova Catalog](https://sne.space) [[GitHub repo](https://github.com/astrocatalogs/supernovae)]
* [The Open TDE Catalog](https://tde.space) [[GitHub repo](https://github.com/astrocatalogs/tidaldisruptions)]
* [The Open Nova Catalog](https://opennova.space) [[GitHub repo](https://github.com/astrocatalogs/novae)]
* [The Open Black Hole Catalog](https://holes.space) [[GitHub repo](https://github.com/astrocatalogs/blackholes)]

The process for creating one's own open astronomy catalog involves checking out this package and designing a "module" for it that is specific to that catalog's needs, a [Wiki is available](https://github.com/astrocatalogs/astrocats/wiki) with instructions for doing so. At the moment the most developed module is the Open Supernova Catalog module; to set up astrocats with the supernovae module, one needs to check out two repositories:

```shell
git clone git@github.com:astrocatalogs/astrocats.git
cd astrocats
git clone git@github.com:astrocatalogs/supernovae.git
```

## Install

This project requires [Python 3.11](https://www.python.org/downloads/) or later and uses [uv](https://docs.astral.sh/uv/) for package management.

```shell
uv python pin 3.11
uv sync
uv run astrocats setup
```

Run the test suite with:

```shell
uv run pytest
```
