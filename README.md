# hcpython
A automated python-based workflow for reproducible and scalable minimal processing with the HCPPipelines.

## Installation

### FSL
```
python fslinstaller.py
```

### Connectome workbench
```bash
mamba install bioconda::connectome-workbench
```

## Running the workflow

```bash
snakemake --cores 8
```