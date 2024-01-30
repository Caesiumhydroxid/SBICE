mutable struct CurrentSource
    I
    terminals::Vector{String}
    name::String
end

function writeToMatrix(A::Matrix, device::CurrentSource, nodeNames::Dict{String,Int}, lastSolution::Vector,  frequency::Float64=0.0) 
end

function RHS(A::Vector, device::CurrentSource, nodeNames::Dict{String,Int}, lastSolution::Vector, frequency::Float64=0.0)
    n1 = nodeNames[device.terminals[1]]
    n2 = nodeNames[device.terminals[2]]
    A[n1] -= device.I
    A[n2] += device.I
end

function amountOfMatrixEntries(device::CurrentSource)
    return 2
end

function nameOfMatrixEntries(device::CurrentSource)
    return device.terminals
end