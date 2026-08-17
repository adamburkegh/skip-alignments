# Skip Probability Computation
In process mining, alignments are a core concept to synchronize actual process executions with a process model.
This repository contains the code to compute _skip probabilities_ for a given log and process tree. We provide all source code and references to replicate the results of skip probabilities from the paper "Skip Probabilities for Subprocesses".

## Library Usage
The project can also be installed and used as a Python library. 
After installing the package, the skip probability functionality can be
imported directly with `from skipalignments import *`, 

## Requirements
- `Python ≥ 3.10`
- `pm4py`
- `pandas`
- `tqdm`
- `Ebi`

## Installation
The project is packaged with `setuptools` (see `pyproject.toml`), with the
package source living under `src/skipalignments`.

```bash
pip install -e .
```

This installs the `skipalignments` package and its dependencies
(`pandas`, `pm4py`, `tqdm`).

## Structure of This Repository
This repository contains everything needed to compute skip alignments and to recreate the evaluation from the paper.

You can recreate the models and probabilities with `im_models.ipynb`, `indulpet_models.ipynb`, and `rand_models.ipynb`. This might take a few days.

The .py files carry the algorithms to compute skip alignments and skip probabilities, and they box the PM4py calls.

## Running Example
An introduction to compute skip probabilities for subprocesses is given in `examples.ipynb`. Is discusses the tree and log from the running example in the paper.

## Required Event Logs
You need to download the event logs used in this repository to recreate the evaluation results. Download, extract, and save the .xes files to disk. You need to provide the paths to these files in each notebook.

- Road Fines: [Download](https://doi.org/10.4121/uuid:270fd440-1057-4fb9-89a9-b699b47990f5)
- Request For Payment: [Download](https://doi.org/10.4121/uuid:895b26fb-6f25-46eb-9e48-0dca26fcd030)
- International Declarations: [Download](https://data.4tu.nl/datasets/91fd1fa8-4df4-4b1a-9a3f-0116c412378f)

## Precomputed Results
We provide the computational results used in our evaluation in the folders `im_results`, `indulpet_results`, and `rand_results`. They are equivalent to the files obtained by running the three notebooks again.

## Ebi
Querying the stochastic path languages in the derivation process requires Ebi. Follow the instructions of [Ebi](https://ebitools.org/) to setup the environment. For skip probabilities, we expect `ebi.exe` to be located in _./src/skipprobabilities_ folder.

## Third Party Dependencies
As scientific library in the Python ecosystem, we rely on external libraries to offer our features. We refer to [this](https://github.com/process-intelligence-solutions/pm4py/tree/release/third_party) page for a detailed list of licenses of the dependencies used in this project. We specifically modified the PM4py library to perform our computations and refer to the PM4py [license](https://github.com/process-intelligence-solutions/pm4py/blob/release/LICENSE) for details.
We additionally refer to [this](https://ebitools.org/) page for licensing information of Ebi.
Finally, we compute skip alignments whose code is yet not publically available. We add reference eonce available.
