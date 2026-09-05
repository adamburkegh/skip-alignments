# Skip Probability Computation
In process mining, alignments are a core concept to synchronize actual process executions with a process model.
This repository contains the code to compute _skip probabilities_ for a given log and process tree. We provide all source code and references to replicate the results of skip probabilities from the paper "Skip Probabilities for Subprocesses".

This implements the techniques from the following two papers:

*Philipp Bär, Sander J.J. Leemans, Moe T. Wynn (2025). A Full Picture in Conformance Checking: Efficiently Summarizing All Optimal Alignments. International Conference on Business Process Management (BPM) 2025.*

*Philipp Bär, Adam T. Burke, Moe T. Wynn, Sander J.J. Leemans (2025). Skip Probabilities for Subprocesses. International Conference on Process Mining (ICPM) 2025.*

## Library Usage
The project can also be installed and used as a Python library. 
After installing the package, the skip alignments and skip 
probability functionality can be imported directly with 
`from skip_alignments import *`, 

## Requirements
- `Python ≥ 3.10`
- `Ebi` - requires separate installation
- Libraries listed under pyproject.toml (immediate dependencies) and 
requirements.txt (output of pip freeze)

## Installation
The project is packaged with `setuptools` (see `pyproject.toml`), with the
package source living under `src/skip_alignments`.

```bash
pip install -e .
```

This installs the `skip-alignments` package and its dependencies
(`pandas`, `pm4py`, `tqdm`).

## Structure of This Repository
This repository contains everything needed to compute skip alignments and to recreate the evaluation from the paper.

You can recreate the models and probabilities with `im_models.ipynb`, `indulpet_models.ipynb`, and `random_models.ipynb`. This might take a few days.

The .py files carry the algorithms to compute skip alignments and skip probabilities, and they box the PM4py calls.

## Running Example
An introduction to compute skip probabilities for subprocesses is given in the notebook `examples/example.ipynb`. Is discusses the tree and log from the running example in the paper. A pure Python version of the same example is in `examples/quickstart.py`.

## Required Event Logs
You need to download the event logs used in this repository to recreate the evaluation results. Download, extract, and save the .xes files to disk. You need to provide the paths to these files in each notebook.

- Road Fines: [Download](https://doi.org/10.4121/uuid:270fd440-1057-4fb9-89a9-b699b47990f5)
- Request For Payment: [Download](https://doi.org/10.4121/uuid:895b26fb-6f25-46eb-9e48-0dca26fcd030)
- International Declarations: [Download](https://data.4tu.nl/datasets/91fd1fa8-4df4-4b1a-9a3f-0116c412378f)

## Precomputed Results
We provide the computational results used in our evaluation in the folders `im_results`, `indulpet_results`, and `rand_results`. They are equivalent to the files obtained by running the three notebooks again.

## Ebi
Querying the stochastic path languages in the derivation process requires Ebi. Follow the instructions of [Ebi](https://ebitools.org/) to setup the environment. For skip probabilities, we expect `ebi` or `ebi.exe` to be available on the 
PATH. The exact path used can be overridden by setting the variable 
`skip_alignments.probabilities.EBI_EXECUTABLE`.

The Ebi CLI's argument shape has changed between versions; the calls in
`probabilities.py` are known to work against Ebi 0.3.14. If Ebi rejects a
call with a usage error, its CLI arguments may have moved on since.


## Third Party Dependencies and licenses
This project is licensed under the AGPL 3.0. This is mainly due to the dependency on pm4py, which is an AGPL project. This project uses pm4py as a library, and includes adapted code from that source.

The code for computing alignments in [alignall.py](src/skip_alignments/alignall.py) is adapted from pm4py, specifically

`pm4py/algo/conformance/alignments/petri_net/variants/state_equation_a_star.py`




