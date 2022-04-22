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
    signs = sign.(Symmetric(get_attribute_layer_weights(attribute, attr)))
    net = generate_complete_network(3)
    triads = get_undir_triads(net)

    params = Params()
    params.attr = attribute
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
end