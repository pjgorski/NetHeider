using DrWatson
quickactivate("NetHeider")
using DataFrames
using NetHeider
using JLD2

foldername = "sims_add_unbalanced"
res = collect_results(datadir(foldername); black_list = ["trans_table", "trans_mean", "trans_std"])

fname = datadir(foldername, "collected_results.jld2")
jldopen(fname, "w") do file
    file["res"] = res
end