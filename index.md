{% include gtag.js %}

## Bioinformatics Toolbox

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
- The size of the image is about **43 GB**.
- Pulling the image (option #1) takes about **15 minutes** on a Google Cloud machine type [**e2-standard-2**](https://cloud.google.com/compute/docs/general-purpose-machines#e2_machine_types) in zone [**us-central1-a**](https://cloud.google.com/compute/docs/regions-zones).
- Building the image (option #2) takes about **210 minutes** on a Google Cloud machine type [**e2-standard-2**](https://cloud.google.com/compute/docs/general-purpose-machines#e2_machine_types) in zone [**us-central1-a**](https://cloud.google.com/compute/docs/regions-zones).

## Citation

If you use `bioinformatics-toolbox` in your work, please cite [10.5281/zenodo.8103969](https://doi.org/10.5281/zenodo.8103969)

---
