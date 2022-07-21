using Dates
using Base: @kwdef


abstract type AbstractParams end

mutable struct Params <: AbstractParams
    "Attribute type"
    attr::AbstractAttributes

    "Number of nodes"
    N::Int

    "Network structure (possible ones: complete, ring, NetSense)"
    net_str::String

    "Parameter of network structure. For ring network it is the number of nearest neighbors. For NetSense dataset it will be the starting semester."
    net_str_param::Int

    "Dynamics parameter rewiring or changing attribute for negative links"
    pr_neg::Float64

    "Dynamics parameter rewiring or changing attribute for positive links"
    pr_pos::Float64

    "Dynamics parameter. Rate of adding a link"
    padd::Float64

    "Dynamics parameter. Probability of closing a triad rather than adding a random link. "
    pclose_triad::Float64

    "Dynamics parameter. ABLTD parameter p"
    pn::Float64

    "Flag if there should be a constant rate of adding links or links should be only considered to be added after finding unbalanced triad"
    const_rate_flag::Bool

    "Function for adding edges"
    add_edges::Function

    "Max number of steps"
    step_max::Int

    "After how many steps a balanced level should calculated"
    measure_balance_every_step::Int

    "Number of repetitions"
    repetitions::Int 

    "Date of starting experiment"
    date::Date 

    "Filename of mean values"
    m_filename::String 

    "After this number of seconds, the state of the simulations will be given"
    inform_after::Int 

    "After this ratio of all reps the state of the simulations will be given"
    inform_every::Float64 

    "After this number of seconds, the results will be saved"
    save_after::Int 

end
export Params

function Params(;
    attr_type = OrderedAttributes,
    attr_params = [5, 0.5, 2], 
    N = 3,
    net_str = "complete",
    pr_neg = 0.1,
    pr_pos = 0.1,
    padd = 0.1,
    pclose_triad = 0.,
    pn = 0.5,
    const_rate_flag = false,
    add_edges = add_single_edge2!,
    step_max = 1000,
    measure_balance_every_step = 10,
    repetitions = 100,
    m_filename = "results.csv",
    inform_after = 60,
    inform_every = 0.25,
    save_after = 60
)
    return Params(
        attr_type(attr_params...),
        N,
        net_str,
        0.,
        pr_neg,
        pr_pos,
        padd,
        pclose_triad,
        pn,
        const_rate_flag,
        add_edges,
        step_max,
        measure_balance_every_step,
        repetitions,
        Date(Dates.now()),
        m_filename,
        inform_after, 
        inform_every, 
        save_after
    )
end

function Base.:(==)(par1::Params, par2::Params)
    fields = [:N, :attr]
    
    for field in fields
        if getfield(par1, field) != getfield(par2, field)
            return false
        end
    end
    return true
end

Base.isequal(par1::Params, par2::Params) = par1 == par2

function Base.hash(par1::Params, h::UInt)
    h = hash(par1.N,h)
    h = hash(par1.attr,h)
    return hash(Params,h)
end

# end