![Bioinformatics Toolbox](toolbox.png)

# Bioinformatics Toolbox

`bioinformatics-toolbox` is a [Docker](https://www.docker.com/) container for common bioinformatics and genomics related tools.

The list of the installed tools and packages can be found [here](Tools.md)

## Installation

### Option 1: Pulling from GitHub Registry (Recommended)

The built image can be downloaded as follows:

```bash
sudo docker run -it ghcr.io/ahmedmoustafa/bioinformatics-toolbox /bin/bash
```

### Option 2: Building from the `Dockerfile`

```bash
git clone https://github.com/ahmedmoustafa/bioinformatics-toolbox.git
```

```bash
cd bioinformatics-toolbox/
```

```bash
sudo docker build -t bioinformatics-toolbox .
```

```bash
sudo docker run -it bioinformatics-toolbox
```

### Notes
- The size of the image is about **60 GB**.
- Pulling the image (option #1) takes about **50 minutes** on a Google Cloud machine type [**e2-medium**](https://cloud.google.com/compute/docs/machine-types) in zone [**us-west2-a**](https://cloud.google.com/compute/docs/regions-zones).
- Building the image (option #2) takes about **six hours** on a Google Cloud machine type [**e2-medium**](https://cloud.google.com/compute/docs/machine-types) in zone [**us-west2-a**](https://cloud.google.com/compute/docs/regions-zones).

## Citation

If you use `bioinformatics-toolbox` in your work, please cite DOI [10.5281/zenodo.4936052](https://doi.org/10.5281/zenodo.4936052)

---
