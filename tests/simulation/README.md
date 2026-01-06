# Simulation instructions

To enable full replication, simulation was run with a Python virtual environment established by `uv` using the process detailed below.
```
uv venv --python 3.13.9
source .venv/bin/activate
uv pip install chromax==0.0.4a0
```

Environment solving led to the following packages being installed:

```
 + chromax==0.0.4a0
 + jax==0.8.2
 + jaxlib==0.8.2
 + jaxtyping==0.3.4
 + ml-dtypes==0.5.4
 + numpy==2.4.0
 + opt-einsum==3.4.0
 + pandas==2.3.3
 + python-dateutil==2.9.0.post0
 + pytz==2025.2
 + scipy==1.16.3
 + six==1.17.0
 + tzdata==2025.3
 + wadler-lindig==0.1.7
```

Matplotlib will subsequently installed with `uv pip install matplotlib` which gave the following package installations:

```
 + contourpy==1.3.3
 + cycler==0.12.1
 + fonttools==4.61.1
 + kiwisolver==1.4.9
 + matplotlib==3.10.8
 + packaging==25.0
 + pillow==12.1.0
 + pyparsing==3.3.1
```

Simulation was then run with `python simulation_pipeline.py -t 14 -o output`. This script automatically sets the default parameters used for the simulation (i.e., `--snpMB 1000 --genomeLength 10000000 --cmMB 3.0 --offspring 10000`) and sets seed values to ensure identical generations where pseudorandomness is involved.
