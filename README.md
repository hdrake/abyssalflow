abyssalflow
==============================
[![Build Status](https://travis-ci.com/hdrake/abyssalflow.svg?branch=master)](https://travis-ci.com/hdrake/abyssalflow)
[![codecov](https://codecov.io/gh/hdrake/abyssalflow/branch/master/graph/badge.svg)](https://codecov.io/gh/hdrake/abyssalflow)
[![License:MIT](https://img.shields.io/badge/License-MIT-lightgray.svg?style=flt-square)](https://opensource.org/licenses/MIT)

Repository for running and analyzing the Planetary Geostrophic Circulation Model

--------

# loading programming environments

To create the python environment `abyssalflow`, I recommend using `acaconda`, which creates a self-contained and self-consistent `python` programming environment with
```bash
conda env create -f environment.yml
```

To activate the Julia environment `abyssalflow`, install Julia v1.2.0 and launch the Julia REPL in the `abyssalflow` project environment in the `/notebooks/` repository with the following commands:
```bash
cd /notebooks/
julia --project=../../abyssalflow
```

# Jupyterlab

Post-processing of PGCM output is done via `jupyter-lab` using `jupyter` notebooks with a Julia v1.2.0 kernel. To ensure that jupyter recognizes the Julia v1.2.0 kernel, I recommend switch to the Julia REPL with `]` and building `IJulia` with the commands
```julia
] build IJulia
```
and then calling `jupyter-lab` within the `abyssalflow` python environment
```bash
[user]$ conda activate abyssalflow
(abyssalflow)[user]$ jupyter-lab
```

--------

<p><small>Project based on the <a target="_blank" href="https://github.com/jbusecke/cookiecutter-science-project">cookiecutter science project template</a>.</small></p>
