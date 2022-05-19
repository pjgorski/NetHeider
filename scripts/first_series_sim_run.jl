using DrWatson
quickactivate("NetHeider")

using NetHeider
using LightGraphs

using Plots

pns = [0:0.1:1...]

th = 0.3
params = Params(; net_str = "NetSense", attr_params = [8, th, 3], repetitions = 10)
net = NetHeider.generate_network_structure(params)

p= (attr=zeros(params.N, params.attr.g),
    signs=zeros(params.N, params.N), signs_old = zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g),
    hlp=zeros(params.N, params.N))

strdicts = @strdict pns
dicts = dict_list(strdicts)

vals = similar(pns)

for (i, dict) in enumerate(dicts)
    pn = let
        @unpack pns = dict
        pns
    end

    params.pn = pn
    
    balanced_table, balanced_mean, balanced_std, last_val, last_std, 
        trans_table, trans_mean, trans_std, bal_unbal_table, bu_mean, bu_std, bal2bal_mean, unbal2bal_mean = 
        performSimulationRepetitions(params; p = p, savefolder = ["data", "testseries"])

    vals[i] = bal2bal_mean[end]
end

plot(pns, vals)

