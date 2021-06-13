![Bioinformatics Toolbox](toolbox.png)

# Bioinformatics Toolbox

`bioinformatics-toolbox` is a [Docker](https://www.docker.com/) container for common bioinformatics and genomics related tools.

The list of the installed tools and packages can be found [here](Tools.md)

## Installation

### Option 1: Pulling from GitHub Registry (Recommended)

The built image can be downloaded as follows:

`sudo docker run -it ghcr.io/ahmedmoustafa/bioinformatics-tools /bin/bash`

### Option 2: Building from the `Dockerfile`

`git clone https://github.com/ahmedmoustafa/bioinformatics-toolbox.git`

`cd bioinformatics-toolbox/`

`sudo docker build -t bioinformatics-toolbox .`

`sudo docker run -it bioinformatics-toolbox`

### Note
The size of the built image is about 16 GB.

## Citation

If you use `bioinformatics-toolbox` in your work, please cite DOI [10.5281/zenodo.4936052](https://doi.org/10.5281/zenodo.4936052)

[![DOI](https://zenodo.org/badge/375832205.svg)](https://zenodo.org/badge/latestdoi/375832205)

---
