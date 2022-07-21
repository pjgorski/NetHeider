using DrWatson
quickactivate("NetHeider")

import Dates
using NetHeider
using LightGraphs

# using Plots

# pns = [0:0.2:1...]
pns = [0.3, 1/3, 0.5, 5/9, 5/7, 0.7]
ths = [0:15...] ./ 16
padds = [0.01:0.04:0.2...]
pclose_triads = [0, 0.35]
# pr_negs = [0:0.1:1...] 
# pr_poss = [0:0.1:1...]
prs = [0:0.05:0.3...]

params = Params(; net_str = "NetSense", attr_params = [8, 0, 3], const_rate_flag = false, repetitions = 10, inform_after = 3600, 
    measure_balance_every_step = 5, step_max = 400)
net = NetHeider.generate_network_structure(params)

p= (attr=zeros(params.N, params.attr.g),
    signs=zeros(params.N, params.N), signs_old = zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g),
    hlp=zeros(params.N, params.N))

strdicts = @strdict pns ths padds pclose_triads prs
dicts = dict_list(strdicts)

# vals = similar(pns)

curtime = time()
for (i, dict) in enumerate(dicts)
    pn, th, padd, pclose_triad, pr = let
        @unpack pns, ths, padds, pclose_triads, prs = dict
        pns, ths, padds, pclose_triads, prs
    end

    params.pn = pn
    params.attr = OrderedAttributes(8, th, 3)
    params.padd = padd
    params.pclose_triad = pclose_triad
    params.pr_pos = pr
    params.pr_neg = pr
    
    balanced_table, balanced_mean, balanced_std, last_val, last_std, 
        trans_table, trans_mean, trans_std, bal_unbal_table, bu_mean, bu_std, bal2bal_mean, unbal2bal_mean, links_num, triads_num = 
        performSimulationRepetitions(params; p = p, savefolder = ["data", "sims_add_unbalanced2"])

    if time() - curtime > params.inform_after
        global curtime = time()

        cur_date = Dates.now()

        print("$cur_date: " * savename(dict))
    end
end



