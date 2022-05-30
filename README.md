###### based on the data from [BlueBrainProject](https://bbp.epfl.ch/nmc-portal) we reconstruct the neocortical microcolumn
To run code, take some steps.
Go to [BlueBrainProject Data](https://bbp.epfl.ch/nmc-portal/downloads.html) and download L1, L23, L4, L5, L6
```
1)unzip L*, and move to folder AllLayers
2)python3 1_prepare_data.py
3)nrnivmodl mechanisms/
4)mpiexec -n 4 nrniv 2_run.py
5)python3 3_plot.py
```