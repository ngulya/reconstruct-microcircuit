###### based on the data from [BlueBrainProject](https://bbp.epfl.ch/nmc-portal) we reconstruct the neocortical microcolumn
To run code, take some steps
```
sh get_data.sh
python3 1_prepare_data.py
nrnivmodl mechanisms/
mpiexec -n 4 nrniv 2_run.py
python3 3_plot.py
```