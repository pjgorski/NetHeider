using DrWatson
quickactivate(@__DIR__)

using NetHeider
using Test

using LinearAlgebra
using SparseArrays
import LightGraphs

function initialize_attr()
    attribute = OrderedAttributes(3, 0.25, 3)
    attr = [1 1 1;
        1 1 3; 
        1 1 2]
    attr = Float64.(attr)
    signs = sign.(Symmetric(get_attribute_layer_weights(attribute, attr)))
    net = generate_complete_network(3)
    triads = get_undir_triads(net)

    params = Params()
    params.attr = attribute
    return params, attr, signs, triads, net
end

function initialize_netsense_attr()
    params = Params(; net_str = "NetSense", attr_params = [8, 0, 3])
    attr = deepcopy(NetHeider.netSenseAttributes)
    signs = sign.(Symmetric(get_attribute_layer_weights(params.attr, attr)))
    signs[signs .== 0.] .= 1
    signs[diagind(signs)] .= 0.

    net = NetHeider.generate_network_structure(params)
    triads = get_undir_triads(net)

    return params, attr, signs, triads, net
end

@testset "Triad update tests" begin
    @testset "removing link" begin
        params, attr, signs, triads, net = initialize_attr()

        #connection will be surely removed
        params.pr_neg = 1.
        params.pr_pos = 1.
        attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)
        @test length(triads) == 0
    end

    @testset "changing attribute in negative link" begin
        params, attr, signs, triads, net = initialize_attr()

        params.pr_pos = 0.
        params.pr_neg = 0.
        params.pn = 1.

        possible_attrs = [copy(attr), copy(attr)]
        possible_attrs[1][1,3] = 2
        possible_attrs[2][2,3] = 2

        obtained_attrs = []
        for _ in 1:100
            _, attr, signs, triads, net = initialize_attr()
            attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)
            if !(attr in obtained_attrs)
                push!(obtained_attrs, copy(attr))
            end
            @test attr in possible_attrs
        end
        for attr in possible_attrs
            @test attr in obtained_attrs
        end
    end

    @testset "changing attribute in positive link" begin
        params, attr, signs, triads, net = initialize_attr()

        params.pr_pos = 0.
        params.pr_neg = 0.
        params.pn = 0.

        possible_attrs = [copy(attr), copy(attr), 
            copy(attr), copy(attr), copy(attr), copy(attr), 
            copy(attr), copy(attr)]
        possible_attrs[1][1,1] = 2
        possible_attrs[2][1,2] = 2
        possible_attrs[3][2,1] = 2
        possible_attrs[4][2,2] = 2
        possible_attrs[5][3,1] = 2
        possible_attrs[6][3,2] = 2
        possible_attrs[7][3,3] = 1
        possible_attrs[8][3,3] = 3

        obtained_attrs = []
        for _ in 1:1000
            _, attr, signs, triads, net = initialize_attr()
            attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)
            if !(attr in obtained_attrs)
                push!(obtained_attrs, copy(attr))
            end
            @test attr in possible_attrs
        end
        for attr in possible_attrs
            @test attr in obtained_attrs
        end
    end
end

@testset "adding edge test" begin
    @testset "adding a single link" begin
        params, attr, signs, triads, net = initialize_attr()

        #connection will be surely removed
        params.pr_neg = 1.
        params.pr_pos = 1.
        attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)
        
        params.padd = 1.
        NetHeider.add_single_edge!(net, params)
        triads = get_undir_triads(net)

        @test length(triads) == 1
    end

    @testset "adding multiple links" begin
        params = Params()
        net = LightGraphs.Graph(3)
        cmpl = generate_complete_network(3)
        params.padd = 1.

        NetHeider.add_edges!(net, params)

        @test net == cmpl
    end

    @testset "adding a triad" begin
        params, attr, signs, triads, net = initialize_netsense_attr()

        net_copy = deepcopy(net)
        params.padd = 1.
        for _ in 1:10
            NetHeider.add_single_edge2!(net_copy, params)
        end
        triads2 = get_undir_triads(net_copy)

        net_copy = deepcopy(net)
        params.padd = 1.
        params.pclose_triad = 1.
        for _ in 1:10
            NetHeider.add_single_edge2!(net_copy, params)
        end
        triads3 = get_undir_triads(net_copy)

        @test length(triads3) > length(triads2)

    end
end

@testset "calculating triad transitions" begin
    @testset "Delta1 -> Delta 0" begin
        params, attr, signs, triads, net = initialize_attr()

        params.pr_pos = 0.
        params.pr_neg = 0.
        params.pn = 1.

        signs_old = copy(signs)

        attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)

        type_trans, bal_unbal = NetHeider.calculate_triad_transitions(signs_old, signs, triads, triads)

        @test sum(abs.(type_trans)) == 1
        @test type_trans[2,1] == 1
        @test sum(abs.(bal_unbal)) == 1
        @test bal_unbal[2,1] == 1
    end

    @testset "Delta1 -> (Delta2,Delta3)" begin
        params, attr, signs, triads, net = initialize_attr()

        params.pr_pos = 0.
        params.pr_neg = 0.
        params.pn = 0.

        signs_old = copy(signs)

        possible_sum_of_signs = [-2, -6]
        was_observed = [false, false]

        for _ in 1:100
            _, attr, signs, triads, net = initialize_attr()
            attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)

            @test sum(signs) in possible_sum_of_signs

            ind = findfirst(sum(signs) .== possible_sum_of_signs)
            if !was_observed[ind]
                was_observed[ind] = true
                type_trans, bal_unbal = NetHeider.calculate_triad_transitions(signs_old, signs, triads, triads)

                if ind == 1
                    @test sum(abs.(type_trans)) == 1
                    @test type_trans[2,3] == 1
                    @test sum(abs.(bal_unbal)) == 1
                    @test bal_unbal[2,1] == 1
                else
                    @test sum(abs.(type_trans)) == 1
                    @test type_trans[2,4] == 1
                    @test sum(abs.(bal_unbal)) == 1
                    @test bal_unbal[2,2] == 1
                end
            end
        end
    end

    @testset "triad vanishes" begin
        params, attr, signs, triads, net = initialize_attr()

        params.pr_pos = 1.
        params.pr_neg = 1.
        params.pn = 0.5

        signs_old = copy(signs)
        triads_old = copy(triads)

        attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)

        type_trans, bal_unbal = NetHeider.calculate_triad_transitions(signs_old, signs, triads_old, triads)

        @test sum(abs.(type_trans[1:4,1:4])) == 0
        @test sum(abs.(bal_unbal[1:2,1:2])) == 0
    end

    @testset "triad appears" begin
        params, attr, signs, triads, net = initialize_attr()

        #connection will be surely removed
        params.pr_neg = 1.
        params.pr_pos = 1.
        attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)
        
        signs_old = copy(signs)
        triads_old = copy(triads)

        params.padd = 1.
        NetHeider.add_single_edge!(net, params)
        triads = get_undir_triads(net)

        #connection will be surely stay
        params.pr_neg = 0.
        params.pr_pos = 0.
        attr, signs, triads, net = NetHeider.single_update(params, attr, signs, triads, net)

        type_trans, bal_unbal = NetHeider.calculate_triad_transitions(signs_old, signs, triads_old, triads)

        @test sum(abs.(type_trans[1:4,1:4])) == 0
        @test sum(abs.(bal_unbal[1:2,1:2])) == 0
    end
end

@testset "simulation repetitions" begin
    params = Params(; N=10)
    # p = (attr=zeros(params.N, params.attr.g), signs=zeros(params.N, params.N), signs_old = zeros(params.N, params.N), new_attr=zeros(params.N, params.attr.g), hlp=zeros(params.N, params.N), adj_mat=Matrix(adjacency_matrix(net, Float64)))

    performSimulationRepetitions(params; savefolder = "test")
end