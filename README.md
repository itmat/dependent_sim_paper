# Generating Correlated Data for Omics Simulation

This repository gives the source code for generating all figures and documents for our manuscript "Generating Correlation Data for Omics Simulation", see our [biorxiv preprint](https://www.biorxiv.org/content/10.1101/2025.01.31.634335v1).

# Dependencies

To run this, you need a system with the following:

 - Python 3.10
 - Apptainer 1.3

Then a virtual environment with the appropriate python packages must be made:

```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Afterwards, the manuscript can be made using Snakemake and apptainer, see the `run_pipeline.sh` file for an example command that was tailored to our compute cluster.
