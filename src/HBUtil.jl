# some functions used to calculate HB in bilayers
# module HBUtil

# using DifferentialEquations
using LinearAlgebra


# checks whether multiplication of weights in all triads is positive.
# It does not have to be close to +-1. Then returns true.
function is_hb(x::Array{Float64,3}, n)
    xs = sign.(Symmetric3d(x, n))
    for l = 1:size(x, 3)
        v = view(xs, :, :, l)
        if sum(v * v .* v) != (n - 1) * n * (n - 2)
            return false
        end
    end
    return true
end
export is_hb

# Case with only one layer
# checks whether multiplication of weights in all triads is positive.
# It does not have to be close to +-1. Then returns true.
function is_hb(x::Array{Float64,2}, n)
    xs = sign.(Symmetric(x))
    if sum(xs * xs .* xs) != (n - 1) * n * (n - 2)
        return false
    end
    return true
end
export is_hb

# Case with only one layer
# checks whether paradise is achieved.
# Checks whether all weights in upper triangle are positive.
# (They don't have to be close to 1).
function is_paradise(x::Array{Float64,2}, n)
    if all(x .>= 0.0)
        return true
    else
        return false
    end
end
export is_paradise

# Case with only one layer
# checks whether hell is achieved.
# Checks whether all weights in upper triangle are negative.
# (They don't have to be close to 1).
function is_hell(x::Array{Float64,2}, n)
    if all(x .<= 0.0)
        return true
    else
        return false
    end
end
export is_hell

# checks whether multiplication of weights in all triads is positive.
# It does not have to be close to +-1. Then returns true.
# It returns three bool values: HB in the whole network, HB in 1st layer and HB in the 2nd layer
function is_hb2(x::Array{Float64,3}, n)
    xs = sign.(Symmetric3d(x, n))
    hb = true
    hb_l = [true, true]
    for l = 1:size(x, 3)
        v = view(xs, :, :, l)
        if sum(v * v .* v) != (n - 1) * n * (n - 2)
            hb = false
            hb_l[l] = false
        end
    end
    (hb, hb_l...)
end
export is_hb2

# Calculates ratio of balanced triads
function get_balanced_ratio(x::Array{Float64,2}, n)
    xs = sign.(Symmetric(x))
    BN = (sum(xs * xs .* xs) + n * (n - 1) * (n - 2)) / 2 / 6 #number of balanced triads
    return BN / n / (n - 1) / (n - 2) * 6
end
export get_balanced_ratio

function get_balanced_ratio(x::Array{Int,2}, n)
    xs = sign.(Symmetric(x))
    BN = (sum(xs * xs .* xs) + n * (n - 1) * (n - 2)) / 2 / 6 #number of balanced triads
    return BN / n / (n - 1) / (n - 2) * 6
end

function get_balanced_ratio_efficient(xs::Matrix{Float64}, triads_count::Int; hlp=zeros(size(xs)))
    mul!(hlp, xs, xs)
    hlp .*= xs
    BN = sum(hlp) / triads_count / 12 + 0.5
    return BN
end
export get_balanced_ratio_efficient

# Calculates balanced ratio in the case of not complete network. 
# xs is a signed network where there are 0s where there are no links. 
# (If that's not the case, one should run `xs .*= adj_mat` beforehand.) 
# adj_mat is an adjacency matrix. 
function get_balanced_ratio_not_complete(xs::Matrix{Float64}, adj_mat::Matrix{Float64}; hlp=zeros(size(xs)))
    mul!(hlp, adj_mat, adj_mat)
    hlp .*= adj_mat
    maxval = sum(hlp) #doing the above on xs can return values from range [-maxval, maxval]

    mul!(hlp, xs, xs)
    hlp .*= xs

    #because of the mentioned above range, one need to transform the output
    BN = (sum(hlp) + maxval) / 2 / maxval
    return BN
end
export get_balanced_ratio_not_complete

# Case with only one layer.
# Returns the counts of triad types (Delta0, Delta1, Delta2, Delta3)
# Weights in upper triangle don't have to be close to 1.
function get_triad_counts(x::Array{Float64,2}, n)
    Deltas = zeros(4)

    xs = sign.(x) .< 0 # gets Boolean array with True, where a link is negative

    for i = 1:1:n, j = (i+1):1:n, k = (j+1):1:n
        # triad = [xs[i,j],xs[i,k],xs[j,k]]

        Deltas[xs[i, j]+xs[i, k]+xs[j, k]+1] += 1
    end

    return Deltas
end
export get_triad_counts

# calculates similarity parameter between the layers
# works only for 2 layers.
function get_layer_similarity(x::Array{Float64,3}, n)
    xs = sign.(triu3d(x, 1))
    s = sum(xs[:, :, 1] .* xs[:, :, 2]) / (n * (n - 1) / 2)
end
export get_layer_similarity

# calculates similarity parameter between two layers (or layer and similarity matrix)
function get_similarity(x1::Array{Float64,2}, x2::Array{Float64,2}, n)
    xs1 = sign.(triu(x1, 1))
    xs2 = sign.(triu(x2, 1))

    s = sum(xs1 .* xs2) / (n * (n - 1) / 2)
end
export get_similarity

# calculates correlation parameter between two layers, where
# xpm1 is the layer where values are +-1 and x2 is the other layer.
# The alg. is my alg.: it is the sum of absolute differences of sign(x2) and xpm1
# weighted by the abs(x2) and normalized by sum(abs(x2)).
function get_correlation(xpm1::Array{Float64,2}, x2::Array{Float64,2})
    u = (triu(xpm1, 1))
    L = (triu(x2, 1))

    s = 1 - sum(abs.(sign.(L) .- sign.(u)) .* abs.(L)) / sum(abs.(L)) / 2 #2 is because sign(L)-u returns 0 or 2
end
export get_correlation

# Returns: 
# * transition Matrix Delta[i,j] indicates the number of times triad of type (i-1) became a triad of type (j-1)
# * transition Matrix of (balanced, unbalanced). 
function calculate_triad_transitions(xs_old::Matrix{Float64}, xs_new::Matrix{Float64}, triads_old, triads_new; hlp=zeros(size(xs))) 
    Deltas = zeros(4,4)
    bal_unbal = zeros(2,2)
    for triad in triads_old
        if !(triad in triads_new)
            continue
        end
        triad_links_inds = ((triad[1], triad[2]), (triad[2], triad[3]), (triad[3], triad[1]))
        triad_links_old = [xs_old[inds...] for inds in triad_links_inds]
        triad_links_new = [xs_new[inds...] for inds in triad_links_inds]

        is_balanced_old_ind = prod(triad_links_old) > 0 ? 1 : 2
        is_balanced_new_ind = prod(triad_links_new) > 0 ? 1 : 2

        bal_unbal[is_balanced_old_ind, is_balanced_new_ind] += 1

        negs_old_ind = sum(triad_links_old .== -1) + 1
        negs_new_ind = sum(triad_links_new .== -1) + 1

        Deltas[negs_old_ind, is_balanced_new_ind] += 1
    end

    return Deltas, bal_unbal
end
# end
