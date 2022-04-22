module NetHeider

using DrWatson

include(srcdir("Attributes.jl"))
include(srcdir("HBUtil.jl"))
include(srcdir("Params.jl"))
include(srcdir("NetworkStructure.jl"))

# plik symulacje
# dostaje na wejście plik params
# generuje zadaną ilość razy atrybuty
# liczy statystyki
# zapisuje do zadanego pliku wyjściowego
include(srcdir("Simulation.jl"))

# w skryptach wywołanie symulacji

end