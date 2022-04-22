using DrWatson

using StatsBase
using LightGraphs

DrWatson.default_prefix(params::Params) = "Experiment_" * string(params.date)
DrWatson.default_allowed(::Params) = (Real, String, AbstractAttributes)
DrWatson.allaccess(::Params) = (:attr, :N, :pr_neg, :pr_pos, :padd, :pn, :net_str)

function single_update(params::Params, attr::Matrix{Int}, signs::Matrix{Float64}, triads, net)
    if length(triads) == 0
        return attr, signs, triads
    end
    #choose triad
    triad = rand(triads)
    triad_links = (signs[triad[1], triad[2]], signs[triad[2], triad[3]], signs[triad[3], triad[1]])
    triad_links_inds = ((triad[1], triad[2]), (triad[2], triad[3]), (triad[3], triad[1]))

    #determine triad type
    is_balanced = prod(triad_links) > 0

    if !is_balanced
        s = sum(triad_links)
        if s == 1 # triad with 1 neg link
            if rand() < params.pn #change neg link
                neg_link_ind = findfirst(triad_links .== -1)
                change_link = triad_links_inds[neg_link_ind]
            else
                pos_link_inds = findall(triad_links .== 1)
                change_link = triad_links_inds[rand(pos_link_inds)]
            end
        else # triad with 3 neg links
            change_link = rand(triad_links_inds)
        end

        pr = signs[change_link...] > 0 ? params.pr_pos : params.pr_neg
        if rand() < pr
            # remove connection
            rem_edge!(net, change_link...)

            triads = get_undir_triads(net)
        else #change attribute
            #choose agent
            ind_1 = rand(1:2)
            agent_1 = change_link[ind_1]
            agent_2 = ind_1 == 1 ? change_link[2] : change_link[1]

            #find set of attributes
            v = get_degeneracy(params.attr)
            dif = attr[agent_1, :] .- attr[agent_2, :]
            if signs[change_link...] > 0 #this link is positive. We want to make it negative
                #find set of not completely different attributes
                attr_inds = findall(abs.(dif) .!= v - 1)

                #exclude attributes that cannot become more different (for instance attribute 1 cannot become smaller)
                attr_inds = [attr_ind for (i, attr_ind) in enumerate(attr_inds) if dif[i] == 0 || (attr[agent_1, attr_ind] != 1 && attr[agent_1, attr_ind] != v)]
            else #this link is negative. We want to make it positive
                #find set of not the same attributes
                attr_inds = findall(attr[agent_1, :] .!= attr[agent_2, :])
            end

            attr_ind = rand(attr_inds)
            if signs[change_link...] > 0 #this link is positive. We want to make it negative

                if dif[attr_ind] == 0
                    if attr[agent_1, attr_ind] == 1
                        possible_new_vals = [attr[agent_1, attr_ind] + 1]
                    elseif attr[agent_1, attr_ind] == v
                        possible_new_vals = [attr[agent_1, attr_ind] - 1]
                    else
                        possible_new_vals = [attr[agent_1, attr_ind] - 1, attr[agent_1, attr_ind] + 1]
                    end
                else
                    possible_new_vals = [attr[agent_1, attr_ind] + sign(dif[attr_ind])]
                end
            else #this link is negative. We want to make it positive
                possible_new_vals = [attr[agent_1, attr_ind] - sign(dif[attr_ind])]
            end

            attr[agent_1, attr_ind] = rand(possible_new_vals)
            signs = sign.(Symmetric(get_attribute_layer_weights(params.attr, attr)))
        end
    end
    return attr, signs, triads, net
end

#add a link with rate padd
#TODO I dont like it, because no new triads will be added
function add_single_edge!(net, params::Params;
    hlp=(adj_mat=Matrix(adjacency_matrix(net, Float64)),))

    if rand() < params.padd
        hlp.adj_mat .= Matrix(adjacency_matrix(net, Float64))

        if sum(1.0 .- hlp.adj_mat) - params.N < 0.5 * params.N^2
            possible_edges = findall(triu(1.0 .- hlp.adj_mat, 1) .== 1)
            chosen_edge = rand(possible_edges).I
        else
            chosen_edge = rand(1:params.N, 2)
            while has_edge(net, chosen_edge...) || chosen_edge[1] == chosen_edge[2]
                chosen_edge = rand(1:params.N, 2)
            end
        end
        add_edge!(net, chosen_edge...)
    end
end
export add_single_edge!

#each unconnected pair is connected with prob padd
#TODO I dont like it, because no new triads will be added
function add_edges!(net, params::Params;
    hlp=(adj_mat=Matrix(adjacency_matrix(net, Float64)),))

    hlp.adj_mat .= Matrix(adjacency_matrix(net, Float64))
    for inds in CartesianIndices(hlp.adj_mat)
        if hlp.adj_mat[inds] == 1
            continue
        elseif inds[1] >= inds[2]
            continue
        end
        if rand() < params.padd
            add_edge!(net, inds.I...)
        end
    end
end
export add_edges!

# p = (attr = zeros(params.N, params.attr.g), weights = zeros(params.N, params.N), signs = zeros(params.N, params.N), new_attr = zeros(params.N, params.attr.g), weights_row = zeros(params.N), signs_row = zeros(params.N), alpha_mod = zeros(params.N), hlp = zeros(params.N, params.N), adj_mat = Matrix(adjacency_matrix(net, Float64)))
# balanced = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)))
function performSimulation!(balanced, params::Params, net=generate_network_structure(params); triads=get_undir_triads(net),
    p=(attr=zeros(params.N, params.attr.g),
        signs=zeros(params.N, params.N), signs_old = copy(signs), new_attr=zeros(params.N, params.attr.g),
        hlp=zeros(params.N, params.N), adj_mat=Matrix(adjacency_matrix(net, Float64)))
)
    attr, signs, new_attr, hlp, adj_mat = p
    #generate network
    # net = generate_network_structure(params)

    #generate attributes
    attr .= get_attributes(params.attr, params.N)
    # weights .= Symmetric(get_attribute_layer_weights(params.attr, attr))

    #generate signed connections (Jij)
    signs .= sign.(Symmetric(get_attribute_layer_weights(params.attr, attr)))

    new_attr .= attr

    #prepare result tables
    # balanced = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)))

    measure_balance_counter = 1
    measure_balance_index = 1

    #run Dynamics
    changes = 0
    for i in 1:params.step_max
        signs_old .= copy(signs)
        triads_old = copy(triads)
        #update triad
        attr, signs, triads, net = single_update(params, attr, signs, triads, net)

        #add a link with rate padd
        params.add_edges(net, params; hlp=p)
        triads = get_undir_triads(net)

        if measure_balance_counter >= params.measure_balance_every_step
            #measure and store the level of balance
            balanced[measure_balance_index] = get_balanced_ratio_not_complete(signs .* adj_mat, adj_mat; hlp=hlp)
            measure_balance_index += 1
            measure_balance_counter = 1
        else
            measure_balance_counter += 1
        end
    end

    if measure_balance_index < size(balanced, 1)
        balanced[measure_balance_index] = get_balanced_ratio_not_complete(signs .* adj_mat, adj_mat; hlp=hlp)
        measure_balance_index += 1
    end

    # print(changes)

    return balanced, changes
    # #save results

    # plp = Array{Float64,1}(undef, params.repetitions)
    # ND = params.N * (params.N - 1) * (params.N - 2) / 6

    # for i = 1:params.repetitions
    #     plp[i] = get_random_polarization(params; attr = attr, weights = weights, triads = triads, ND = ND)
    # end

    # #save Results
    # @tagsave(
    #     datadir("sims", savename(params, "jld2")),
    #     @strdict params plp
    # )
end
export performSimulation!

function performSimulationRepetitions(params::Params; p=(attr=zeros(params.N, params.attr.g),
        weights=zeros(params.N, params.N), signs=zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g),
        weights_row=zeros(params.N), signs_row=zeros(params.N), alpha_mod=zeros(params.N), hlp=zeros(params.N, params.N)),
    savefolder="")

    attr, weights, signs, new_attr, weights_row, signs_row, alpha_mod, hlp = p

    #generate network
    #Assuming the network topology is not random!
    net = generate_network_structure(params)

    p = (p..., adj_mat=Matrix(adjacency_matrix(net, Float64)))

    #prepare result tables
    balanced_table = zeros(params.repetitions, Int(ceil(params.step_max / params.measure_balance_every_step)))
    balanced_mean = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)))
    balanced_std = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)))

    for rep in 1:params.repetitions
        #run dynamics
        bal_row = @view balanced_table[rep, :]
        performSimulation!(bal_row, params, net; p=p)
    end

    balanced_mean .= mean(balanced_table, dims=1)[:]
    balanced_std .= std(balanced_table, dims=1)[:]

    last_val = balanced_mean[end]
    last_std = balanced_std[end]

    #saving Results
    @tagsave(
        datadir("sims", savefolder, savename(params, "jld2")),
        @strdict params balanced_table balanced_mean balanced_std last_val last_std
    )

    return balanced_table, balanced_mean, balanced_std, last_val, last_std
end
export performSimulationRepetitions

function get_random_polarization(params::Params; attr=zeros(params.N, params.attr.g), weights=zeros(params.N, params.N), triads=zeros(4), ND=params.N * (params.N - 1) * (params.N - 2) / 6)
    attr .= get_attributes(params.attr, params.N)
    weights .= get_attribute_layer_weights(params.attr, attr)

    triads .= get_triad_counts(weights, params.N)

    return (triads[3] + triads[4]) / ND
end

#calculates Hamiltonian given a signed array `signs` and similarity weights `weights`. 
#It calculates for node of index `index` or if `index==0` for the whole network. 
function calc_hamiltonian(params::Params, signs, weights, index)
    if index == 0
        return sum([calc_hamiltonian(params, signs, weights, i) for i in 1:params.N])
    end

    alpha_mod = [Jij * (Jij > 0 ? params.alpha : 1.0 - params.alpha) for Jij in signs[index, :]]
    Hi = sum(-alpha_mod .* weights[index, :])

    return Hi
end

#calculates Hamiltonian given a signed row (for a given node) `signs_row` and similarity weights `weights_row`.
function calc_hamiltonian(params::Params, signs_row, weights_row; p=(alpha_mod=zeros(params.N),))

    p.alpha_mod .= ifelse.(signs_row .> 0, -params.alpha, signs_row .* (params.alpha - 1.0))
    p.alpha_mod .*= weights_row
    Hi = sum(p.alpha_mod)

    return Hi
end