using Printf
using TOML
using HDF5
using Dumper
using CairoMakie

@doc raw"""
prune_analysis(data_file::AbstractString, observable_tags::Vector{String}, collapse_beta::Float64)

The MC data of a given observable is plotted for different seeds and user enters start and end;
The final data is written in a HDF5 File containing everything;
As an input needs to provide the path of the file where the data is stored. 
We need to provide the tag of the observables to be visualized and of the weights.
"""
function prune_analysis(data_file::AbstractString, observable_tags::Vector{String}, collapse_beta::Float64)
    # Open HDF5 file
    h5file = h5open(data_file, "r")

    # Read betas and find closest index
    beta_list = read(h5file["betas_traj"])[:]
    target_index = argmin(abs.(beta_list .- collapse_beta))
    @show target_index

    # Plot observables
    f = Figure()
    for (i, tag) in enumerate(observable_tags)
        data = read(h5file[tag])  # assume shape (trajectory, beta/time)
        vals = data[target_index, :]

        ax = Axis(f[i, 1], title=tag, xlabel="Time Series", ylabel=tag)
        scatter!(ax, 1:length(vals), vals)
        lines!(ax, 1:length(vals), vals)
    end
    #display(f)
    save("tmp.png", f)
    
    # Ask once for start/end
    print("Global Start Index (default 1): ")
    rd = readline()
    rd_stripped = strip(rd)
    start = rd_stripped != "" ? parse(Int, rd_stripped) : 1

    print("Global End Index (default end): ")
    rd = readline()
    rd_stripped = strip(rd)
    # Assume all observables have same time dimension length change this in the future, if we do this we need to be carefull with the weights we take!
    len = size(read(h5file[observable_tags[1]]), 2)
    @show len
    endd = rd_stripped != "" ? parse(Int, rd_stripped) : len

    start = max(start, 1)
    endd = min(endd, len)
    
    println("Saving data to new file")
    # Save pruned data to a new HDF5 file
    dir = dirname(data_file)
    pruned_file = joinpath(dir,"outfile.pruned.hdf5")
    
    outfile = DumpFile(pruned_file)

    # assure that we write different seeds into the same final file!
    if !isfile(pruned_file)
        dump!(outfile,"betas_traj",beta_list)
    end

    # Copy weights pruned!

    @show size(read(h5file["norm_metts"]))
    
    weights = read(h5file["norm_metts"])[:,start:endd]
    
    dump!(outfile, "norm_metts", weights)
    
    for tag in observable_tags
        data = read(h5file[tag])
        pruned = data[:, start:endd]
        dump!(outfile, tag, pruned)
    end
    close(h5file)

    println("Pruned data written to:")
    println(pruned_file)
end
