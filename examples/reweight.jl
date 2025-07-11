using Mbar4METTS

path = "triangularJ1J2.L.12.W.3.J2.0.500/Beta_Collapse.1.000/tau.0.20.cutoff.1.0e-06.maxdim.512/outfile.pruned.hdf5"

# Final array to compute the temperatures
beta_final = LinRange(0,10,200)


# Get_energie
reweight_observable(beta_final,path,"energy")

# Get_dE/DB
reweight_first_dev_temp_observable(beta_final,path,"energy")