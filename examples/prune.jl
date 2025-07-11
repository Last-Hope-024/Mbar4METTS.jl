using Mbar4METTS

path = "triangularJ1J2.L.12.W.3.J2.0.500/Beta_Collapse.1.000/tau.0.20.cutoff.1.0e-06.maxdim.512/outfile.seed.1.hdf5"
prune_analysis(path,["energy","entanglement"],1.0)