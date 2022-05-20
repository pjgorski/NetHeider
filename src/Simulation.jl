using DrWatson

using StatsBase
using LightGraphs

DrWatson.default_prefix(params::Params) = "Experiment_" * string(params.date)
DrWatson.default_allowed(::Params) = (Real, String, AbstractAttributes)
DrWatson.allaccess(::Params) = (:attr, :N, :pr_neg, :pr_pos, :padd, :pn, :pclose_triad, :net_str)

function single_update(params::Params, attr::Matrix{Float64}, signs::Matrix{Float64}, triads, net)
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
            old_attr = deepcopy(attr)
            #choose agent
            ind_1 = rand(1:2)
            agent_1 = change_link[ind_1]
            agent_2 = ind_1 == 1 ? change_link[2] : change_link[1]

            #find set of attributes
            v = get_degeneracy(params.attr)
            dif = attr[agent_1, :] .- attr[agent_2, :]
            if signs[change_link...] > 0 #this link is positive. We want to make it negative
                #find set of not completely different attributes
                # attr_inds = findall(abs.(dif) .!= v - 1)

                # #exclude attributes that cannot become more different (for instance attribute 1 cannot become smaller)
                # attr_inds = [attr_ind for (i, attr_ind) in enumerate(attr_inds) if dif[i] == 0 || (attr[agent_1, attr_ind] != 1 && attr[agent_1, attr_ind] != v)]
                attr_inds = [i for i in 1:length(dif) if dif[i] == 0 || (attr[agent_1, i] != 1 && attr[agent_1, i] != v)]
            else #this link is negative. We want to make it positive
                #find set of not the same attributes
                attr_inds = findall(attr[agent_1, :] .!= attr[agent_2, :])
            end

            if length(attr_inds) == 0
                # display(attr)
                # display(triad)
                # display(v)
                # display(dif)
                # display((agent_1, agent_2))
            else
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
                signs[signs .== 0.] .= 1
                signs[diagind(signs)] .= 0.
            end

            if !are_attributes_correct(params.attr, attr)
                display(attr)
                display(triad)
                display(v)
                display(dif)
                display((agent_1, agent_2))
                display(attr_inds)
                display(attr_ind)
                display(possible_new_vals)
                display(attr[[agent_1, agent_2], :])
                display(old_attr[[agent_1, agent_2], :])
                display(net)
                display(length(triads))
                error("wrong attributes!")
            end
        end
    end
    return attr, signs, triads, net, is_balanced
end
export single_update

#add a link with rate padd
#TODO I dont like it, because no new triads will be added
function add_single_edge!(net, params::Params;
    hlp=(adj_mat=Matrix(adjacency_matrix(net, Float64)),))

    if rand() < params.padd
        hlp.adj_mat .= Matrix(adjacency_matrix(net, Float64))

        if sum(1.0 .- hlp.adj_mat) - params.N < 0.5 * params.N^2
            possible_edges = findall(triu(1.0 .- hlp.adj_mat, 1) .== 1)
            if length(possible_edges) == 0
                return 
            end
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

# adds one link with prob padd. With prob pclose_triad this link closes a triad. 
function add_single_edge2!(net, params::Params;
        hlp=(adj_mat=Matrix(adjacency_matrix(net, Float64)),))

    if rand() < params.padd
        hlp.adj_mat .= Matrix(adjacency_matrix(net, Float64))
        if rand() < params.pclose_triad
            wedges_t, _ = get_wedges_single_edges(hlp.adj_mat)
            wedge = rand(wedges_t)
            if hlp.adj_mat[wedge[1], wedge[2]] == 0
                chosen_edge = wedge[1], wedge[2]
            elseif hlp.adj_mat[wedge[1], wedge[3]] == 0
                chosen_edge = wedge[1], wedge[3]
            else chosen_edge = wedge[2], wedge[3]
            end
        else
            if sum(1.0 .- hlp.adj_mat) - params.N < 0.5 * params.N^2
                possible_edges = findall(triu(1.0 .- hlp.adj_mat, 1) .== 1)
                if length(possible_edges) == 0
                    return 
                end
                chosen_edge = rand(possible_edges).I
            else
                chosen_edge = rand(1:params.N, 2)
                while has_edge(net, chosen_edge...) || chosen_edge[1] == chosen_edge[2]
                    chosen_edge = rand(1:params.N, 2)
                end
            end
            add_edge!(net, chosen_edge...)

            hlp.adj_mat .= Matrix(adjacency_matrix(net, Float64))
        end
    end
end
export add_single_edge2!

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

# p = (attr=zeros(params.N, params.attr.g), signs=zeros(params.N, params.N), signs_init = zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g), hlp=zeros(params.N, params.N), adj_mat=Matrix(adjacency_matrix(net, Float64)))
# res = (balanced = zeros(Int(ceil(params.step_max / params.measure_balance_every_step))), triad_trans = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)), 4, 4), bal_unbal = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)), 2, 2))
function performSimulation!(res, params::Params, net=generate_network_structure(params); triads = [],
    p=(attr=zeros(params.N, params.attr.g),
        signs=zeros(params.N, params.N), signs_init = zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g),
        hlp=zeros(params.N, params.N), adj_mat = zeros(params.N, params.N))
)
    if isempty(triads) 
        triads=get_undir_triads(net)
    end
    if isempty(p.adj_mat)
        p.adj_mat .= Matrix(adjacency_matrix(net, Float64))
    end
    attr, signs, signs_init, new_attr, hlp, adj_mat = p
    
    #generate attributes
    if params.net_str == "NetSense"
        attr .= NetHeider.netSenseAttributes
    else
        attr .= get_attributes(params.attr, params.N)
    end
    # weights .= Symmetric(get_attribute_layer_weights(params.attr, attr))

    #generate signed connections (Jij)
    signs .= sign.(Symmetric(get_attribute_layer_weights(params.attr, attr)))
    signs[signs .== 0.] .= 1
    signs[diagind(signs)] .= 0.

    new_attr .= attr

    #prepare result tables
    # balanced = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)))
    # signs_old .= copy(signs)
    # triads_old = copy(triads)
    triads_init = copy(triads)
    signs_init .= copy(signs)

    measure_balance_counter = 1
    measure_balance_index = 1

    #run Dynamics
    for i in 1:params.step_max
        #update triad
        # print("u")
        # print(net.ne)
        attr, signs, triads, net, was_triad_balanced = single_update(params, attr, signs, triads, net)
        # print("a")
        if params.const_rate_flag || !was_triad_balanced
            #add a link with rate padd
            params.add_edges(net, params; hlp=p)
            triads = get_undir_triads(net)
        end

        if measure_balance_counter >= params.measure_balance_every_step
            #measure and store the level of balance
            adj_mat .= Matrix(adjacency_matrix(net, Float64))
            res.balanced[measure_balance_index] = get_balanced_ratio_not_complete(signs .* adj_mat, adj_mat; hlp=hlp)

            # (balanced = bal_row, signs_row = signs_row, triads_row = triads_row, links_row = links_row, triads_num = triads_num)

            if !isempty(res.signs_row)
                res.signs_row[measure_balance_index, :, :] .= signs
            end
            if !isempty(res.triads_row)
                res.triads_row[measure_balance_index] = triads
            end
            res.links_row[measure_balance_index] = net.ne
            res.triads_num[measure_balance_index] = length(triads)

            res.triad_trans[measure_balance_index, :, :], 
                res.bal_unbal[measure_balance_index, :, :] = calculate_triad_transitions(signs_init, signs, triads_init, triads)

            # res.triad_trans[measure_balance_index, :, :], 
            #     res.bal_unbal[measure_balance_index, :, :] = calculate_triad_transitions(signs_old, signs, triads_old, triads)
            # signs_old .= copy(signs)
            # triads_old = copy(triads)

            measure_balance_index += 1
            measure_balance_counter = 1
        else
            measure_balance_counter += 1
        end
    end

    if measure_balance_index < size(res.balanced, 1)
        res.balanced[measure_balance_index] = get_balanced_ratio_not_complete(signs .* adj_mat, adj_mat; hlp=hlp)
        # res.triad_trans[measure_balance_index, :, :], 
        #     res.bal_unbal[measure_balance_index, :, :] = calculate_triad_transitions(signs_old, signs, triads_old, triads)

        # res.signs_row[measure_balance_index, :, :] .= signs
        # res.triads_row[measure_balance_index] = triads
        if !isempty(res.signs_row)
            res.signs_row[measure_balance_index, :, :] .= signs
        end
        if !isempty(res.triads_row)
            res.triads_row[measure_balance_index] = triads
        end
        res.links_row[measure_balance_index] = net.ne
        res.triads_num[measure_balance_index] = length(triads)

        res.triad_trans[measure_balance_index, :, :], 
            res.bal_unbal[measure_balance_index, :, :] = calculate_triad_transitions(signs_init, signs, triads_init, triads)

        measure_balance_index += 1
    end

    # print(changes)

    return res
end
export performSimulation!

function performSimulationRepetitions(params::Params; p=(attr=zeros(params.N, params.attr.g),
    signs=zeros(params.N, params.N), signs_old = zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g),
    hlp=zeros(params.N, params.N)),
    savefolder=["data","sims"])

    attr, signs, signs_old, new_attr, hlp = p

    #generate network
    #Assuming the network topology is not random!
    net = generate_network_structure(params)
    p = (p..., adj_mat=Matrix(adjacency_matrix(net, Float64)))

    #prepare result tables
    balanced_table = zeros(params.repetitions, Int(ceil(params.step_max / params.measure_balance_every_step)))
    balanced_mean = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)))
    balanced_std = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)))

    trans_table = zeros(params.repetitions, Int(ceil(params.step_max / params.measure_balance_every_step)), 4, 4)
    bal_unbal_table = zeros(params.repetitions, Int(ceil(params.step_max / params.measure_balance_every_step)), 2, 2)
    trans_mean = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)), 4, 4)
    trans_std = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)), 4, 4)
    bu_mean = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)), 2, 2)
    bu_std = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)), 2, 2)

    # signs_table = zeros(Int, params.repetitions, Int(ceil(params.step_max / params.measure_balance_every_step)), params.N, params.N)
    # triads_table = Array{Array{Any, 1}, 2}(undef, params.repetitions, Int(ceil(params.step_max / params.measure_balance_every_step)))

    links_num = zeros(Int, params.repetitions, Int(ceil(params.step_max / params.measure_balance_every_step)))
    triads_num = zeros(Int, params.repetitions, Int(ceil(params.step_max / params.measure_balance_every_step)))

    for rep in 1:params.repetitions
        # display("Started " * string(rep))
        net = generate_network_structure(params)
        p.adj_mat .= Matrix(adjacency_matrix(net, Float64))

        #run dynamics
        bal_row = @view balanced_table[rep, :]
        trans_row = view(trans_table, rep, :, :, :)
        bal_unbal_row = view(bal_unbal_table, rep, :, :, :)
        # signs_row = view(signs_table, rep, :, :, :)
        # triads_row = view(triads_table, rep, :, :, :)
        signs_row = []
        triads_row = []
        links_row = view(links_num, rep, :)
        triads_num_row = view(triads_num, rep, :)

        performSimulation!((balanced = bal_row, triad_trans = trans_row, bal_unbal = bal_unbal_row, 
            signs_row = signs_row, triads_row = triads_row, links_row = links_row, triads_num = triads_num_row), 
            params, net; p=p)
    end

    balanced_mean .= mean(balanced_table, dims=1)[:]
    balanced_std .= std(balanced_table, dims=1)[:]

    trans_mean .= mean(trans_table, dims=1)[1,:,:,:]
    trans_std .= std(trans_table, dims=1)[1,:,:,:]

    bu_mean .= mean(bal_unbal_table, dims=1)[1,:,:,:]
    bu_std .= std(bal_unbal_table, dims=1)[1,:,:,:]

    bal2bal_mean = bu_mean[:,1,1] ./ sum(bu_mean[:,1,:], dims=2)
    unbal2bal_mean = bu_mean[:,2,1] ./ sum(bu_mean[:,2,:], dims=2)

    last_val = balanced_mean[end]
    last_std = balanced_std[end]

    threshold = params.attr.threshold
    what_to_save = @strdict params threshold links_num triads_num balanced_table balanced_mean balanced_std last_val last_std trans_table trans_mean trans_std bal_unbal_table bu_mean bu_std bal2bal_mean unbal2bal_mean
    for field in fieldnames(typeof(params))
        val = getfield(params, field)
        what_to_save[String(field)] = val
    end

    if !isa(savefolder, Array)
        savefolder = [savefolder]
    end

    #saving Results
    safesave(
        projectdir(savefolder..., savename(params, "jld2")),
        what_to_save
    )

    return balanced_table, balanced_mean, balanced_std, last_val, last_std, trans_table, trans_mean, trans_std, bal_unbal_table, bu_mean, bu_std, bal2bal_mean, unbal2bal_mean
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