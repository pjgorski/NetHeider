using DrWatson
quickactivate("NetHeider")
using DataFrames
using NetHeider
using JLD2

foldername = "sims_add_unbalanced2"
res = collect_results(datadir(foldername); black_list = ["trans_table", "trans_mean", "trans_std"], rinclude = [r"Experiment*"])

fname = datadir(foldername, "collected_results2.jld2")
jldopen(fname, "w") do file
    file["res"] = res
end

include(projectdir("scripts", "notebook_res_anal_funs.jl"))
res = load(datadir(foldername, "collected_results2.jld2"))["res"]

pn = [0.4, 0.3, 1/3, 0.5, 5/9, 5/7, 0.7, 1., 0.6]
dict = @dict pn
part_res = get_part_dataframe(res, dict; verbose = true)

fname = datadir(foldername, "collected_partial_results.jld2")
jldopen(fname, "w") do file
    file["res"] = res
end