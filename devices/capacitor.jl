mutable struct Capacitor
    C
    terminals::Vector{String}
    name::String
end

function writeToMatrix(A::Matrix, device::Capacitor, nodeNames::Dict{String,Int},  lastSolution::Vector, frequency::Float64=0.0) 
    n1 = nodeNames[device.terminals[1]]
    n2 = nodeNames[device.terminals[2]]
    A[n1, n1] += device.C * frequency *2 * pi * 1im
    A[n1, n2] -= device.C * frequency *2 * pi * 1im
    A[n2, n1] -= device.C * frequency *2 * pi * 1im
    A[n2, n2] += device.C * frequency *2 * pi * 1im
end

function RHS(A::Vector, device::Capacitor, nodeNames::Dict{String,Int}, lastSolution::Vector, frequency::Float64=0.0)
    #Do nothing
end

function writeToSymMatrix(A::Matrix, device::Capacitor, nodeNames::Dict{String,Int},  lastSolution::Vector, frequency::Float64=0.0) 
    n1 = nodeNames[device.terminals[1]]
    n2 = nodeNames[device.terminals[2]]
    A[n1, n1] += device.name*".C*frequency*1im "
    A[n1, n2] -= device.name*".C*frequency*1im "
    A[n2, n1] -= device.name*".C*frequency*1im "
    A[n2, n2] += device.name*".C*frequency*1im "
end

function nameOfMatrixEntries(device::Capacitor)
    return device.terminals
end