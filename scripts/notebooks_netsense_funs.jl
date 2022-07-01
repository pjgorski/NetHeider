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
using PyCall

py"""
import pickle
 
def load_pickle(fpath):
    with open(fpath, "rb") as f:
        data = pickle.load(f)
    return data
"""

load_pickle = py"load_pickle"

py"""
import networkx as nx
def get_triads(net):
    A = nx.adjacency_matrix(net).todense() 

    N = A.shape[0]

    triads = []

    for i in range(0,N):
        for j in range(i+1,N):
            if A[i,j] > 0:
                for k in range(j+1,N):
                    if A[i,k] > 0 and A[j,k] > 0:
                        triads.append((i,j,k))
    
    return triads

def get_adjacency_matrix(net, ordering = None):
    return nx.adjacency_matrix(net, ordering).todense() 
"""

function get_triads2(A)
    

    N = size(A)[1]

    triads = []

    for i in 1:N
        for j in (i+1):N
            if A[i,j] > 0
                for k in (j+1):N
                    if A[i,k] > 0 && A[j,k] > 0
                        push!(triads, (i,j,k))
                    end
                end
            end
        end
    end
    
    return triads
end

get_triads = py"get_triads"
get_adjacency_matrix = py"get_adjacency_matrix"

function get_all_nodes(file)
    nodes = [n for n in file[1].nodes]
    for i in 2:6
        nodes = [nodes..., [n for n in file[i].nodes()]...]
    end
    return unique(nodes)
end

function create_series_adj_mat(file, nodes = get_all_nodes(file))
    As = [get_adjacency_matrix(file[i], nodes) for i in 1:6]
    return As
end

function get_wedges_single_edges(A)
    N = size(A)[1]

    wedges = []
    single_edges = []

    for i in 1:N
        for j in (i+1):N
            for k in (j+1):N
                if A[i,k] + A[j,k] + A[i,j] == 2
                    push!(wedges, (i,j,k))
                elseif A[i,k] + A[j,k] + A[i,j] == 1
                    push!(single_edges, (i,j,k))
                end
            end
        end
    end
    
    return wedges, single_edges
end

function understand_triad_apperance(As)
    num_mats = length(As)

    all_triads = [get_triads2(As[i]) for i in 1:num_mats]
    all_wedges = []
    all_singles = []
    for i in 1:num_mats
        wedges, single_edges = get_wedges_single_edges(As[i])
        push!(all_wedges, wedges)
        push!(all_singles, single_edges)
    end

    

    triad_appear = zeros(num_mats,num_mats)
    triad_appear_wedges = zeros(num_mats,num_mats)
    triad_appear_single = zeros(num_mats,num_mats)
    triad_appear_nowhere = zeros(num_mats,num_mats)

    for i in 1:num_mats
        for j in i:num_mats
            triad_appear[i,j] = length(all_triads[j]) - length(intersect(all_triads[i], all_triads[j]))
            new_triads = setdiff(all_triads[j], all_triads[i])
            @assert length(new_triads) == triad_appear[i,j]
            triad_appear_wedges[i,j] = length(intersect(new_triads, all_wedges[i]))

            triad_appear_single[i,j] = length(intersect(new_triads, all_singles[i]))
            # triad_vanish_abs[i,j] = length(all_triads[i]) - length(intersect(all_triads[i], all_triads[j]))
            triad_appear_nowhere[i,j] = triad_appear[i,j] - triad_appear_wedges[i, j] - triad_appear_single[i, j]
        end
    end
    return triad_appear, triad_appear_wedges, triad_appear_single, triad_appear_nowhere
end

# defining functions to calibrate this probability
function add_edges!(A1r, edges_to_add, close_triad)
    N = size(A1r)[1]

    added = 0
    while added < edges_to_add
        wedges_t, single_t = get_wedges_single_edges(A1r)
        if rand() < close_triad
            wedge = rand(wedges_t)
            if A1r[wedge[1], wedge[2]] == 0
                i, j = wedge[1], wedge[2]
            elseif A1r[wedge[1], wedge[3]] == 0
                i, j = wedge[1], wedge[3]
            else i, j = wedge[2], wedge[3]
            end
        else
            i, j = rand(1:N, 2)
            if i==j continue
            elseif A1r[i,j] > 0
                continue
            end
        end
        A1r[i,j] = A1r[j,i] = 1
        added += 1
    end
end

function calc_difs(close_triad_prob; term_pairs = [(1, 2), (2, 3), (2, 4), (3, 4), (4, 5), (5, 6)], repetitions = 10)
    dif = 0

    for _ in 1:repetitions
        for term_pair in term_pairs
            start_term = term_pair[1]
            inc = term_pair[2] - start_term

            A1r = As2[start_term] .* As2[start_term + inc] 
            add_edges!(A1r, appear_abs_triads[start_term, start_term + inc], close_triad_prob)
            triad_appear3, triad_appear_wedges3, triad_appear_single3, triad_appear_nowhere3 = understand_triad_apperance([As2[start_term], A1r])

            dif += abs(triad_appear[term_pair...] - triad_appear3[1,2])
            dif += abs(triad_appear_wedges[term_pair...] - triad_appear_wedges3[1,2])
            dif += abs(triad_appear_single[term_pair...] - triad_appear_single3[1,2])
            dif += abs(triad_appear_nowhere[term_pair...] - triad_appear_nowhere3[1,2])
        end
    end
    return dif
end

function calc_difs_only_numT(close_triad_prob; term_pairs = [(1, 2), (2, 3), (2, 4), (3, 4), (4, 5), (5, 6)], repetitions = 10)
    dif = 0

    for _ in 1:repetitions
        for term_pair in term_pairs
            start_term = term_pair[1]
            inc = term_pair[2] - start_term

            A1r = As2[start_term] .* As2[start_term + inc] 
            add_edges!(A1r, appear_abs_triads[start_term, start_term + inc], close_triad_prob)
            triad_appear3, triad_appear_wedges3, triad_appear_single3, triad_appear_nowhere3 = understand_triad_apperance([As2[start_term], A1r])

            dif += abs(triad_appear[term_pair...] - triad_appear3[1,2])
        end
    end
    return dif
end

# a change is a change of attribute from 1 to 2. a change from 1 to 3 are two changes. 
function count_attr_changes(ord_attributes)
    sems = length(ord_attributes)
    N = size(ord_attributes[1])[1]

    changes = zeros(sems, sems)
    max_change = zeros(sems, sems)

    # for sem1 in 1:sems, sem2 in 1:sems
    for i in 1:sems
        for j in i:sems
            attr_change = abs.(ord_attributes[i] .- ord_attributes[j])

            changes[i,j] = sum(filter(!isnan, attr_change))
            max_change[i,j] = length(filter(!isnan, attr_change))*2
        end
    end
    return changes, max_change
end

# a change is a change of attribute from 1 to 2. a change from 1 to 3 are two changes. 
function count_attr_stays(ord_attributes)
    sems = length(ord_attributes)
    N = size(ord_attributes[1])[1]

    stays = zeros(sems, sems)
    max_stays = zeros(sems, sems)

    # for sem1 in 1:sems, sem2 in 1:sems
    for i in 1:sems
        for j in i:sems
            attr_stay = sum(ord_attributes[i] .== ord_attributes[j])

            stays[i,j] = attr_stay
            max_stays[i,j] = length(filter(!isnan, ord_attributes[i] .+ ord_attributes[j]))
        end
    end
    return stays, max_stays
end

#returns two values: 
# first value (true- balanced, false - unbalanced)
# second value is triad type (0, 1, 2, 3) - number of negative links
function get_triad_type(triad, signs)
    links = [(triad[1], triad[2]), (triad[3], triad[2]), (triad[1], triad[3])]

    triad_signs = [signs[link...] for link in links]

    is_balanced = prod(triad_signs) > 0
    triad_type = sum(triad_signs .== -1)

    return is_balanced, triad_type
end

# Returns two Arrays. 
# First - balanced (1), unbalanced (2) or lack of triad (3). 
# Second - triad_type (0-3) or lack of triad (4). 
function get_triads_history(all_triads_list, all_signs; sem_count = 6, all_triads = all_triads)
    T = length(all_triads_list)

    triads_balanced_hist = Array{Any}(undef, T, sem_count)
    triads_type_hist = Array{Any}(undef, T, sem_count)
    for (ti, triad) in enumerate(all_triads_list), sem in 1:sem_count
        is_balanced = 3
        triad_type = 4
        if triad in all_triads[sem]
            is_balanced, triad_type = get_triad_type(triad, all_signs[sem])
            is_balanced = is_balanced ? 1 : 2
        end
        triads_balanced_hist[ti, sem] = is_balanced
        triads_type_hist[ti, sem] = triad_type
    end

    return triads_balanced_hist, triads_type_hist
end

function get_triads_history(threshold::Float64, attrs, all_triads; sem_count = 6)
    aT = OrderedAttributes(8, threshold, 3)

    all_triads_list = unique([triad for triads in all_triads for triad in triads])
    
    all_signs = [Symmetric(sign.(get_attribute_layer_weights(aT, attrs[i]))) .* 1. for i in 1:6]
    map(x-> x[x .== 0.] .= 1, all_signs)

    return get_triads_history(all_triads_list, all_signs; sem_count = sem_count, all_triads = all_triads)
end


function get_triad_transitions(all_triads, threshold, attrs; all_triads_list = [] )
    aT = OrderedAttributes(8, threshold, 3)

    if isempty(all_triads_list)
        all_triads_list = unique([triad for triads in all_triads for triad in triads])
    end
    T = length(all_triads_list)

    all_signs = [Symmetric(sign.(get_attribute_layer_weights(aT, attrs[i]))) .* 1. for i in 1:6]
    map(x-> x[x .== 0.] .= 1, all_signs)

    triads_balanced_hist, triads_type_hist = get_triads_history(all_triads_list, all_signs)

    bal_trans = zeros(3,3) 
    bal_trans2 = zeros(2,2) # let's remove all nans, so if the triad exist in semesters 1, 3, 4, then there are 2 transitions 1->3, 3->4

    type_trans = zeros(5,5)
    type_trans2 = zeros(4,4)

    for (ti, triad) in enumerate(all_triads_list)
        # for sem in 1:5
        #     inds = triads_balanced_hist[ti, sem:(sem+1)]
        #     bal_trans[inds...] += 1
        #     inds = triads_type_hist[ti,sem:(sem+1)] .+ 1
        #     type_trans[inds...] += 1
        # end
        inds = [[triads_balanced_hist[ti,i], triads_balanced_hist[ti,i+1]] for i in 1:5]
        map((x) -> bal_trans[x...] += 1, inds)

        inds = [[triads_type_hist[ti,i], triads_type_hist[ti, i+1]] .+ 1 for i in 1:5]
        map((x) -> type_trans[x...] += 1, inds)

        triad_exists = triads_balanced_hist[ti,:] .!= 3
        els = sum(triad_exists)
        if els == 1
            continue
        end
        
        bal_hist = triads_balanced_hist[ti, triad_exists]
        type_hist = triads_type_hist[ti,triad_exists]

        inds = [[bal_hist[i], bal_hist[i+1]] for i in 1:els-1]
        map((x) -> bal_trans2[x...] += 1, inds)

        inds = [[type_hist[i], type_hist[i+1]] .+ 1 for i in 1:els-1]
        map((x) -> type_trans2[x...] += 1, inds)

    end

    return bal_trans, type_trans, bal_trans2, type_trans2
end

function calc_delta1_vanishing_triad_fate(all_triads, threshold, attrs; all_triads_list = [] )
    aT = OrderedAttributes(8, threshold, 3)

    if isempty(all_triads_list)
        all_triads_list = unique([triad for triads in all_triads for triad in triads])
    end
    T = length(all_triads_list)

    all_signs = [Symmetric(sign.(get_attribute_layer_weights(aT, attrs[i]))) .* 1. for i in 1:6]
    map(x-> x[x .== 0.] .= 1, all_signs)

    _, triads_type_hist = get_triads_history(all_triads_list, all_signs; all_triads = all_triads)

    pos2not = 0
    neg2not = 0
    # cnt = 0
    for sem in 1:5
        for (ind, triad) in enumerate(all_triads_list)
            if triads_type_hist[ind, sem] == 1
                if triads_type_hist[ind, sem + 1] == 4 #it vanished
                    if any(isnan.(ord_attributes[sem+1][collect(triad),:]))
                        continue
                    else
                        i,j,k = triad
                        links = [[i,j], [i,k], [j,k]]

                        vanished = [As[sem+1][link...] == 0 for link in links]

                        # display(triad)
                        # display(links)
                        # display([As[sem+1][link...] for link in links])
                        # display(vanished)
                        
                        #we're intersted only in triads that becaome wedges
                        if sum(vanished) > 1
                            continue
                        end

                        l_ind = findfirst(vanished)

                        sign_changes = [all_signs[sem+1][link...] != all_signs[sem][link...] for link in links]
                        sign_changes[l_ind] = 0

                        # we're interested only in triads that one edge vanished. 
                        if sum(sign_changes) > 0 
                            continue
                        end

                        if all_signs[sem][links[l_ind]...] == 1
                            pos2not += 1
                        else
                            neg2not += 1
                        end

                        # cnt += 1
                        # display((ind, sem))
                    end
                end
            end
        end
    end

    return pos2not, neg2not
end

function get_all_links2(As)
    N = size(As[1])[1]

    all_links = sum(As) .> 0
    links = []

    for i in 1:N
        for j in (i+1):N
            if all_links[i,j] > 0
                push!(links, (i,j))
            end
        end
    end
    
    return links
end

function get_links_history(all_links_list, all_signs)
    D = length(all_links_list)

    links_sign_hist = Array{Any}(undef, D, 6)
    for (ti, link) in enumerate(all_links_list), sem in 1:6
        is_pos = 3
        l_sign = all_signs[sem][link...]
        if isnan(l_sign)
            is_pos = 3
        elseif l_sign == 1
            is_pos = 1
        else
            is_pos = 2
        end
        links_sign_hist[ti, sem] = is_pos
    end

    return links_sign_hist
end


function get_link_transitions(As, threshold, attrs; all_links_list = [] )
    aT = OrderedAttributes(8, threshold, 3)

    if isempty(all_links_list)
        all_links_list = get_all_links2(As)
    end
    D = length(all_links_list)

    all_signs = [Symmetric(sign.(get_attribute_layer_weights(aT, attrs[i]))) .* 1. for i in 1:6]
    map(x-> x[x .== 0.] .= 1, all_signs)

    links_sign_hist = get_links_history(all_links_list, all_signs)

    sign_trans = zeros(3,3) 
    sign_trans2 = zeros(2,2) # let's remove all nans, so if the triad exist in semesters 1, 3, 4, then there are 2 transitions 1->3, 3->4

    for (ti, link) in enumerate(all_links_list)
        # for sem in 1:5
        #     inds = triads_balanced_hist[ti, sem:(sem+1)]
        #     bal_trans[inds...] += 1
        #     inds = triads_type_hist[ti,sem:(sem+1)] .+ 1
        #     type_trans[inds...] += 1
        # end
        inds = [[links_sign_hist[ti,i], links_sign_hist[ti,i+1]] for i in 1:5]
        map((x) -> sign_trans[x...] += 1, inds)

        link_exists = links_sign_hist[ti,:] .!= 3
        els = sum(link_exists)
        if els == 1
            continue
        end
        
        sign_hist = links_sign_hist[ti, link_exists]

        inds = [[sign_hist[i], sign_hist[i+1]] for i in 1:els-1]
        map((x) -> sign_trans2[x...] += 1, inds)

    end

    return sign_trans, sign_trans2
end


function get_sign_trans_probs(all_links, threshold, attrs; all_links_list = [])
    res = get_link_transitions(all_links, threshold, attrs; all_links_list = all_links_list)

    pos_flip = res[1][1,1] / sum(res[1][1,1:2])
    neg_flip = res[1][2,1] / sum(res[1][2,1:2])
    pos_flip2 = res[2][1,1] / sum(res[2][1,1:2])
    neg_flip2 = res[2][2,1] / sum(res[2][2,1:2])

    return pos_flip, neg_flip, pos_flip2, neg_flip2
end

function edge_randomization(all_links, threshold, attrs = ord_attributes, starting_adj_term = 1; pos_to_pos_prob = NaN, neg_to_pos_prob = NaN, all_links_list = [])
    aT = OrderedAttributes(8, threshold, 3)
    
    all_signs = [Symmetric(sign.(get_attribute_layer_weights(aT, attrs[i]))) .* 1. for i in 1:6]
    map(x-> x[x .== 0.] .= 1, all_signs)

    if isnan(pos_to_pos_prob) || isnan(neg_to_pos_prob)
        pos_to_pos_prob, neg_to_pos_prob, pos_flip2, neg_flip2 = get_sign_trans_probs(all_links, threshold, attrs; all_links_list = all_links_list)
    end
    
    starting_adj = As[starting_adj_term]
    starting_signs = all_signs[starting_adj_term]

    N = size(starting_adj)[1]
    new_signs = zeros(N, N)
    for i in 1:N
        for j in i+1:N
            if starting_adj[i, j] > 0
                if starting_signs[i, j] == 1
                    is_pos = rand() < pos_to_pos_prob
                else
                    is_pos = rand() < neg_to_pos_prob
                end
                new_signs[i, j] = 2*is_pos - 1
                new_signs[j, i] = new_signs[i, j]
            end
        end
    end
    return new_signs, starting_signs, pos_to_pos_prob, neg_to_pos_prob
end

function get_triad_transitions2(triads, all_signs )
    sem_count = length(all_signs)

    T = length(triads)

    triads_balanced_hist, triads_type_hist = get_triads_history(triads, all_signs; sem_count = sem_count, all_triads = repeat([triads], sem_count))

    bal_trans = zeros(3,3) 
    bal_trans2 = zeros(2,2) # let's remove all nans, so if the triad exist in semesters 1, 3, 4, then there are 2 transitions 1->3, 3->4

    type_trans = zeros(5,5)
    type_trans2 = zeros(4,4)

    for (ti, triad) in enumerate(triads)
        # for sem in 1:5
        #     inds = triads_balanced_hist[ti, sem:(sem+1)]
        #     bal_trans[inds...] += 1
        #     inds = triads_type_hist[ti,sem:(sem+1)] .+ 1
        #     type_trans[inds...] += 1
        # end
        inds = [[triads_balanced_hist[ti,i], triads_balanced_hist[ti,i+1]] for i in 1:sem_count-1]
        map((x) -> bal_trans[x...] += 1, inds)

        inds = [[triads_type_hist[ti,i], triads_type_hist[ti, i+1]] .+ 1 for i in 1:sem_count-1]
        map((x) -> type_trans[x...] += 1, inds)

        triad_exists = triads_balanced_hist[ti,:] .!= 3
        els = sum(triad_exists)
        if els == 1
            continue
        end
        
        bal_hist = triads_balanced_hist[ti, triad_exists]
        type_hist = triads_type_hist[ti,triad_exists]

        inds = [[bal_hist[i], bal_hist[i+1]] for i in 1:els-1]
        map((x) -> bal_trans2[x...] += 1, inds)

        inds = [[type_hist[i], type_hist[i+1]] .+ 1 for i in 1:els-1]
        map((x) -> type_trans2[x...] += 1, inds)

    end

    return bal_trans, type_trans, bal_trans2, type_trans2
end

function get_attr_transitions(attrs )
    sem_count = length(attrs)

    N, G = size(attrs[1])

    attr_change_count_any = zeros(G) # a change from 1 to 3 is counted as 1
    attr_change_count_abs = zeros(G) # a change from 1 to 3 is counted as 2
    attr_change_count_any2 = zeros(G) # NaNs are discarded
    attr_change_count_abs2 = zeros(G) # NaNs are discarded

    # attr_change_probs_any = attr_change_count_any ./ attr_change_max_any
    # max values of above
    attr_change_max_any = zeros(G) # a change from 1 to 3 is counted as 1
    attr_change_max_abs = zeros(G) # a change from 1 to 3 is counted as 2
    attr_change_max_any2 = zeros(G) # NaNs are discarded
    attr_change_max_abs2 = zeros(G) # NaNs are discarded
    for g = 1:G
        for n in 1:N
            attr_history = [attrs[sem][n,g] for sem in 1:sem_count]

            type_of_change = [attr_history[i+1] - attr_history[i] for i in 1:sem_count-1]

            attr_exists = .!(isnan.(type_of_change))
            els = sum(attr_exists)
            
            attr_change_count_any[g] += sum((type_of_change .!= 0) .& attr_exists)
            attr_change_max_any[g] += els

            attr_change_count_abs[g] += sum(abs.(type_of_change[attr_exists]))
            attr_change_max_abs[g] += 2*els

            attr_exists = .!(isnan.(attr_history))
            els = sum(attr_exists)
            if els == 1
                continue
            end

            attr_no_nan = attr_history[attr_exists]
            type_of_change = [attr_no_nan[i+1] - attr_no_nan[i] for i in 1:els-1]

            attr_change_count_any2[g] += sum(type_of_change .!= 0)
            attr_change_max_any2[g] += els

            attr_change_count_abs2[g] += sum(abs.(type_of_change))
            attr_change_max_abs2[g] += 2*els
        end
    end

    attr_change_probs_any = attr_change_count_any ./ attr_change_max_any
    attr_change_probs_abs = attr_change_count_abs ./ attr_change_max_abs
    attr_change_probs_any2 = attr_change_count_any2 ./ attr_change_max_any2
    attr_change_probs_abs2 = attr_change_count_abs2 ./ attr_change_max_abs2

    return attr_change_probs_any, attr_change_probs_abs, attr_change_probs_any2, attr_change_probs_abs2
end

function attr_randomization(attrs = ord_attributes; attr_change_probs = [], num_of_turns = 1)
    
    if isempty(attr_change_probs) 
        attr_change_probs_any, attr_change_probs_abs, attr_change_probs_any2, attr_change_probs_abs2 = get_attr_transitions(attrs)

        attr_change_probs = attr_change_probs_abs
    end

    probs = attr_change_probs ./ num_of_turns
    
    sem_count = length(attrs)

    N, G = size(attrs[1])

    new_attrs = deepcopy(attrs)
    rand_vals = zeros(G)

    for _ in 1:num_of_turns
        for sem in 1:sem_count, n in 1:N
            if isnan(new_attrs[sem][n,1])
                continue
            end
            rand!(rand_vals)

            change_attr_flag = rand_vals .< probs

            for g in 1:G
                if change_attr_flag[g]
                    if mod(new_attrs[sem][n,g],2) == 1 
                        new_attrs[sem][n,g] = 2
                    elseif rand() < 0.5
                        new_attrs[sem][n,g] = 1
                    else
                        new_attrs[sem][n,g] = 3
                    end
                end
            end
        end
    end
    return new_attrs
end


# a change is a change of attribute from 1 to 2. a change from 1 to 3 are two changes. 
function count_attr_changes_cor_sign(attrs, all_triads)
    sems = length(attrs)
    N = size(attrs[1])[1]

    changes_all_triads = []

    # for sem1 in 1:sems, sem2 in 1:sems
    for i in 1:sems-1
        j = i+1

        changes_triads = zeros(length(all_triads[i]))

        attr_change = abs.(attrs[i] .- attrs[j])
        for (ti, triad) in enumerate(all_triads[i])

            if any(isnan.(attr_change[[triad...],:])) #missing attributes
                changes_triads[ti] = NaN
            else
                changes_triads[ti] = sum(attr_change[[triad...], :])
            end
        end

        push!(changes_all_triads, changes_triads)
    end
    return changes_all_triads
end

function get_balanced_ratio(all_triads, threshold, attrs; all_triads_list = [] )

    aT = OrderedAttributes(8, threshold, 3)

    if isempty(all_triads_list)
        all_triads_list = unique([triad for triads in all_triads for triad in triads])
    end
    T = length(all_triads_list)

    all_signs = [Symmetric(sign.(get_attribute_layer_weights(aT, attrs[i]))) .* 1. for i in 1:6]
    map(x-> x[x .== 0.] .= 1, all_signs)

    triads_balanced_hist, _ = get_triads_history(all_triads_list, all_signs; all_triads = all_triads)

    bal_ratio = [sum(triads_balanced_hist[:,i] .== 1) / (T - sum(triads_balanced_hist[:,i] .== 3)) for i in 1:6]

    all_bal_ratio = sum(triads_balanced_hist .== 1) / (6*T - sum(triads_balanced_hist .== 3))

    return bal_ratio, all_bal_ratio
end

function form_triad_network(all_triads, triad_nodes; all_triads_list = [])
    if isempty(all_triads_list)
        all_triads_list = unique([triad for triads in all_triads for triad in triads])
    end

    N2 = length(triad_nodes)
    N = maximum(triad_nodes)

    As3 = [zeros(N,N) for i in 1:6]

    for sem in 1:6
        for triad in all_triads[sem]
            i, j, k = triad
            As3[sem][i,j] = As3[sem][j,i] = 1
            As3[sem][i,k] = As3[sem][k,i] = 1
            As3[sem][k,j] = As3[sem][j,k] = 1
        end
    end

    return As3
end
