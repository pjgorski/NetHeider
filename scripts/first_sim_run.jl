using DrWatson
quickactivate("NetHeider")

using NetHeider
using LightGraphs

using Plots

th = 0.3
params = Params(; net_str = "NetSense", attr_params = [8, th, 3])
net = NetHeider.generate_network_structure(params)

p= (attr=zeros(params.N, params.attr.g),
    signs=zeros(params.N, params.N), signs_old = zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g),
    hlp=zeros(params.N, params.N))

p = (p..., adj_mat=Matrix(adjacency_matrix(net, Float64)))


bal_row = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)))
trans_row = zeros(Int(ceil(params.step_max / params.measure_balance_every_step)), 4, 4)
bal_unbal_row = zeros( Int(ceil(params.step_max / params.measure_balance_every_step)), 2, 2)

signs_row = zeros(Int, Int(ceil(params.step_max / params.measure_balance_every_step)), params.N, params.N)
triads_row = Array{Array{Any, 1}, 1}(undef, Int(ceil(params.step_max / params.measure_balance_every_step)))

links_row = zeros(Int, Int(ceil(params.step_max / params.measure_balance_every_step)))
triads_num_row = zeros(Int, Int(ceil(params.step_max / params.measure_balance_every_step)))


res = (balanced = bal_row, triad_trans = trans_row, bal_unbal = bal_unbal_row, 
    signs_row = signs_row, triads_row = triads_row, links_row = links_row, triads_num = triads_num_row)

performSimulation!(res, params, net; p = p);

plot(bal_row)
# print(changes)

unbal2bal_mean = res.bal_unbal[:,2,1] ./ sum(res.bal_unbal[:,2,:], dims=2)
bal2bal_mean = res.bal_unbal[:,1,1] ./ sum(res.bal_unbal[:,1,:], dims=2)

plot(bal2bal_mean)
plot!(unbal2bal_mean)