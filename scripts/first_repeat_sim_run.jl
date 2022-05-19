using DrWatson
quickactivate("NetHeider")

using NetHeider
using LightGraphs

# using Plots
th = 0.3
params = Params(; net_str = "NetSense", attr_params = [8, th, 3])
net = NetHeider.generate_network_structure(params)

p= (attr=zeros(params.N, params.attr.g),
    signs=zeros(params.N, params.N), signs_old = zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g),
    hlp=zeros(params.N, params.N))

# p = (p..., adj_mat=Matrix(adjacency_matrix(net, Float64)))

balanced_table, balanced_mean, balanced_std, last_val, last_std = performSimulationRepetitions(params; p = p, savefolder = ["test"])

