# Reproducing pyCSEP: A Software Toolkit for Earthquake Forecast Developers

[![Software](https://zenodo.org/badge/DOI/10.5281/zenodo.6626265.svg)](https://doi.org/10.5281/zenodo.6626265)


A reproducibility package contains the code and data needed to replicate the figures in pyCSEP: An Enhanced Python Toolkit for Earthquake Forecast Developers

We provide the user with options to run Jupyter Notebooks


## Table of contents

* [Instructions](#instructions)
   * [Run the computational environment](#run-the-computational-environment)
      * [Easy-mode using Docker](#easy-mode-using-docker)
      * [Using a conda environment](#using-a-conda-environment)
   * [Reproduce figures individually](#reproduce-figures-individually)
* [Code description](#code-description)
* [Using the Notebooks](#notebook-description)
* [Software versions](#software-versions)
* [Computational effort](#computational-effort)
* [References](#references)


## Instructions
If you have obtained the software from Zenodo, you may skip this step. Make sure that the file contents are extracted and navigate to the package directory.

If you are viewing this on GitHub or have not downloaded the code from Zenodo, open a terminal and download the reproducibility package from GitHub with
```
git clone https://github.com/KennyGraham1/pyCSEP_followup_paper.git
```

Navigate to the newly downloaded directory
```
cd pyCSEP_followup_paper
```

### Run the computational environment

> Important: Use the Docker instructions to ensure the same computational environment as the publication. The `conda` instructions are useful for users that wish to extend the package or want to work with pyCSEP outside of this context.

Now you have two options how to run the package:
 * [Easy-mode using Docker](#easy-mode-using-docker)
 * [Using a conda environment](#using-a-conda-environment)

The easiest way to run the reproducibility package is to run the version of the package in an environment provided by Docker. If you are interested in working with pyCSEP in more detail, we recommend that you install pyCSEP in a `conda` environment.

We have accompanying scripts that work under Linux/macOS.

#### Easy mode using Docker

You will need to have the Docker runtime environment installed and running on your machine. Some instructions
can be found [here](https://docs.docker.com/engine/install/). The following commands will not work unless the Docker engine
is correctly installed on your machine.

If on Linux/maxOS, call in the Terminal/Console:
```
./run_all.sh
```


This step does the following things: (1) download and verify the checksum of the downloaded
data; (2) build a docker image with the computational environment; and (3) run the reproducibility package.

To reproduce the _full_ version, call:
> ```
> ./run_all.sh --full
> ```


When finished, the figures will be stored in `./output`. These can be compared against the expected results that are found in the `expected_results` directory.

The `start_docker.sh` or `start_docker.bat` scripts provide an interactive terminal to re-create individual figures. See [below](#reproduce-figures-individually) for instructions.


> Note: If you are running Docker on MacOS you might run into permission errors when trying to start (or run) the Docker container. To fix this, manually create the `output` (eg `cd output` and add these to the Docker host.)

#### Using a conda environment

Installation instructions can be found in the [pyCSEP documentation](https://docs.cseptesting.org/getting_started/installing.html). **Warning**: this method does not guarantee a stable computing environment, because `conda` provides the newest packages that are compatible with your environment. Please report any issues related to this [here](https://github.com/SCECcode/pycsep/issues).

Create and activate a new conda environment:
```
conda env create -n pycsep_srl
conda activate pycsep_srl
```

Install v0.6.3 of pyCSEP:
```
conda install --channel conda-forge pycsep=0.6.3
```

Download data from Zenodo:
```
./download_data.sh
```

> Note: to download the _'full'_ version, append ` --full` to the command (see [above](#easy-mode-using-docker))

Run the package to reproduce all figures from the manuscript that are supported by your downloaded version
```
cd scripts
python plot_all.py
```

Once completed, the figures can be found in the `outpot` directory in the top-level directory. These can be compared against the expected results that are found in the `expected_results` directory.

To recreate individual figures, follow the instructions below.


## Reproduce figures individually

The scripts to reproduce the figures are contained in the `scripts` directory; change to it if you are in the top-level directory:
```
cd scripts
```
> Note: Any script must be launched from within this `scripts` directory.

Here is an example to recreate Fig. 2:
```
python plot_fig_2.py 
```


## Code description

The top-level directory contains a few helpful scripts for working with the Docker environment. Descriptions of the files in the top-level directory are as follows (the `.sh` and `.bat` scripts provide the same functionality on different operating systems):

* `download_data.{sh|bat}`: downloads and verifies checksums of the data from Zenodo
* `build_docker.{sh|bat}`: (re)builds the Docker image for this environment
* `start_docker.{sh|bat}`: starts Docker container and provides command-line interface with pycsep environment active
* `run_docker.{sh|bat}`: runs the Docker container and automatically launches the package
* `entrypoint.sh`: entrypoint for the runnable Docker container

The code to execute the main experiment can be found in the `scripts` directory of this repository. The files are named
according to the figure they create in the manuscript. The script `plot_all.py` will generate all of the figures supported by the downloaded version.
Descriptions of the files in the `scripts` directory are as follows:

* `plot_fig_1.py`: Plot the Catalog figure
* `plot_fig_2.py`: Plot GEAR1 model projected to New Zealand
* `plot_fig_4to6.py`: Plot figure 4 to 6
            4a.  Pseudo-prospective negative-binomial number test results for four time-invariant seismicity models in New Zealand
            4b. Binary S-test results for earthquake forecasting models in New Zealand during the 2014-2022 evaluation period
            5   Quantitative comparisons of earthquake forecasting models for New Zealand based  on a) Poisson and binary joint log-likelihood scores, b) Kagan I1 information scores and c) Brier scores. 
            6.  a) Receiver Operating Characteristic (ROC) curves, obtained using the alarm-based approach, comparing a Poisson uniform (SUP) seismicity model anddifferent spatially specific time-invariant earthquake forecasting models for New Zealand. 
               b) ROC concentration curves of the models compared with both the SUP model and GeoNet’s catalogue of observed earthquakes 
               c) Comparison of Molchan diagrams depicting the predictive performance of models for  NewZealand, including comparisons among themselves and with SUP.

## Using the Notebooks

To run the notebooks 

```
cd notebooks/
```
> Note: Any script must be launched from within this `scripts` directory.

Here is an example to open and run the jupyter notebook:
```
jupyter-notebook
```
Select of the notbooks and play with it.

## Software versions
* `python>=3.11`
* `pycsep=0.6.3`

To obtain the environment used for publishing this manuscript use [Docker](#easy-mode-using-docker). Advanced users can recreate the environment using `conda` running on Ubuntu 20.04 LTS.


## References

Krafczyk, M. S., Shi, A., Bhaskar, A., Marinov, D., and Stodden, V. (2021).
Learning from reproducing computational results: introducing three principles and the reproduction package.
_Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences, 379_(2197).
doi: [10.1098/rsta.2020.0069](https://doi.org/10.1098/rsta.2020.0069)
