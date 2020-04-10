###### based on the data from [BlueBrainProject](https://bbp.epfl.ch/nmc-portal) we reconstruct the neocortical microcolumn
To run code, take some steps
```
nrniv 1_prepare_data.py
nrnivmodl mechanisms/
nrniv 2_check_data.py
mpiexec -n 3 nrniv 3_run.py
python3 4_plot.py
```