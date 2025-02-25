# Data for Replication

This directory contains the data necessary for replicating the experiments.

The **Sioux Falls** dataset is sourced from the paper *Equilibrium decomposed optimization: a heuristic for the continuous equilibrium network design problem（1987）*

The remaining network datasets are obtained from the [TNTP github repository](https://github.com/bstabler/TransportationNetworks).

Please note that the TNTP repository also includes a version of the **Sioux Falls** dataset. While it appears similar to the version used in the paper after scaling, there are subtle differences between the two datasets. You can use the functions `load_SiouxFalls()` and `load_SiouxFalls2()` to compare them.