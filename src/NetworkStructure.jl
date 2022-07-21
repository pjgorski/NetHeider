using SparseArrays
import LightGraphs
using PyCall

function generate_complete_network(N::Int)
    return LightGraphs.complete_graph(N)
end
export generate_complete_network

function generate_ring_network(N::Int, K::Int)
    return LightGraphs.watts_strogatz(N, K, 0)
end
export generate_ring_network

const netSenseFile = begin
    py"""
    import pickle
    
    def load_pickle(fpath):
        with open(fpath, "rb") as f:
            data = pickle.load(f)
        return data
    """

    load_pickle = py"load_pickle"

    fname = datadir("exp_pro", "NetSense_network.pkl")
    file = load_pickle(fname)

    file
end

const netSenseNetworks, netSenseAttributes_s = begin
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

    def get_adjacency_matrix(net):
        return nx.adjacency_matrix(net).todense() 
    """

    get_triads = py"get_triads"
    get_adjacency_matrix = py"get_adjacency_matrix"

    As = [get_adjacency_matrix(netSenseFile[i]) for i in 1:6]
    net = [LightGraphs.Graph(As[i]) for i in 1:6]

    #creating list of attributes
    list_of_attributes = Set([k for n in netSenseFile[1].nodes for k in keys(netSenseFile[1].nodes[n+1])])
    ord_list_of_attributes = setdiff(list_of_attributes, ["ethnicity_1", "gender_1", "hometown_1", "age_1"])

    #creating table of attribute values
    ord_attributes = Dict((i, zeros(length(netSenseFile[i].nodes), length(ord_list_of_attributes))) for i in 1:6)

    for sem in 1:6
        for (i, node) in enumerate(netSenseFile[sem].nodes())
            for (j, attribute) in enumerate(ord_list_of_attributes)
                ord_attributes[sem][i,j] = netSenseFile[sem].nodes[node + 1][attribute]
            end
        end
    end

    net, ord_attributes
end

const netSenseNetwork, netSenseAttributes = netSenseNetworks[1], netSenseAttributes_s[1]

function generate_NetSenseNetwork(; sem = 1)
    return copy(netSenseNetworks[sem])
end

function generate_network_structure(params::Params)
    if params.net_str == "complete"
        return generate_complete_network(params.N)
    elseif params.net_str == "ring"
        return generate_ring_network(params.N, params.net_str_param)
    elseif params.net_str == "NetSense"
        
        if isempty(params.net_str_param) || params.net_str_param == 0
            sem = 1
        else sem = params.net_str_param
        end
        net = generate_NetSenseNetwork(; sem = sem)

        params.N = length(net.fadjlist)
        params.attr = OrderedAttributes(8, params.attr.threshold, 3)
        
        return net
    end
end
export generate_network_structure

# Returns triads of given network. 
function get_undir_triads(net)
    N = length(net.fadjlist)

    triads = []
    for i in 1:N
        neigh_i = LightGraphs.neighbors(net, i)
        for j in neigh_i
            if j < i
                continue
            end
            neigh_j = LightGraphs.neighbors(net, j)

            ks = intersect(neigh_i, neigh_j)
            # remove smaller indices
            ks = [k for k in ks if k > j]
            for k in ks
                push!(triads, (i, j, k))
            end
        end
    end
    return triads
end
export get_undir_triads

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
export get_wedges_single_edges