using DrWatson
quickactivate(@__DIR__)

using NetHeider
using LinearAlgebra
using LightGraphs
using Plots
using DataFrames
using JLD2
using Statistics
using StatsPlots


# dict should be create the following way. df has fields :a, :b, then:
# a = [4] # where 4 is the value which extracted rows should have
# b = [2]
# dict = @dict a b
function get_part_dataframe(df::DataFrame, dict::Dict; verbose = true)
    cropped_res = deepcopy(df)
    for (field, values) in dict
        inds = findall(x->x in values, cropped_res[:, field])
        cropped_res = cropped_res[inds,:]
    end
    if verbose
        display("Extracted " * string(nrow(cropped_res)) * " rows.")
    end
    return cropped_res
end


function calc_error(b2b1, u2b1, b2b2, u2b2; method = "mse")

    b2b_non_nan = (.!(isnan.(b2b1)) .+ .!(isnan.(b2b2)) ) .== 2

    # if sum(b2b_non_nan) < length(b2b1)
    #     return calc_error(b2b1[b2b_non_nan], u2b1, b2b2[b2b_non_nan], u2b2; method = method)
    # end

    u2b_non_nan = (.!(isnan.(u2b1)) .+ .!(isnan.(u2b2)) ) .== 2

    # if sum(u2b_non_nan) < length(b2b1)
    #     return calc_error(b2b1, u2b1[u2b_non_nan], b2b2, u2b2[u2b_non_nan]; method = method)
    # end

    if method == "mse"
        error = sum((b2b1[b2b_non_nan] .- b2b2[b2b_non_nan]).^2) + sum((u2b1[u2b_non_nan] .- u2b2[u2b_non_nan]).^2)
    elseif method == "mae"
        error = sum(abs.(b2b1[b2b_non_nan] .- b2b2[b2b_non_nan])) + sum(abs.(u2b1[u2b_non_nan] .- u2b2[u2b_non_nan]))
    elseif method == "mse_notbeg"
        # b2b_non_nan[1:3] = 0
        return calc_error(b2b1[4:end], u2b1[4:end], b2b2[4:end], u2b2[4:end]; method = "mse")
    elseif method == "mae_notbeg"
        return calc_error(b2b1[4:end], u2b1[4:end], b2b2[4:end], u2b2[4:end]; method = "mae")
    end

    return error
end

function plot_dict(dicts, df, param_ind)

    dict = dicts[results_df[param_ind, :dict_ind]]
    cropped_res = get_part_dataframe(res, dict)

    time_ind = Int(results_df[param_ind, :time_ind])
    xvals = cropped_res[:, :threshold]
    b2b_vals = [val[time_ind] for val in cropped_res[:, :bal2bal_mean]]
    u2b_vals = [val[time_ind] for val in cropped_res[:, :unbal2bal_mean]]

    p1 = plot(xvals, b2b_vals, markershape = markers[1], label = "ABM T(+->+)")#, linestyle = :none)
    plot!(p1, xvals, u2b_vals, markershape = markers[2], label = "ABM T(-->+)")#, linestyle = :none)
    plot!(p1, real_trans_df.threshold, real_trans_df.netsense_b2b, markershape = markers[3], label = "Net T(+->+)")#, linestyle = :none)
    plot!(p1, real_trans_df.threshold, real_trans_df.netsense_u2b, markershape = markers[4], label = "Net T(-->+)")#, linestyle = :none)
    title!(p1, "steps="*string(10*time_ind))
    
    return p1
end