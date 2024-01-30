mutable struct Resistor
    R
    terminals::Vector{String}
    name::String
end

function writeToMatrix(A::Matrix, device::Resistor, nodeNames::Dict{String,Int},  lastSolution::Vector, frequency::Float64=0.0) 
    n1 = nodeNames[device.terminals[1]]
    n2 = nodeNames[device.terminals[2]]
    A[n1, n1] += 1/device.R
    A[n1, n2] -= 1/device.R
    A[n2, n1] -= 1/device.R
    A[n2, n2] += 1/device.R
end

function RHS(A::Vector, device::Resistor, nodeNames::Dict{String,Int}, lastSolution::Vector, frequency::Float64=0.0)
    #Do nothing
end

function writeToMatrixSym(A::Matrix, device::Resistor, nodeNames::Dict{String,Int},  lastSolution::Vector, frequency::Float64=0.0) 
    n1 = nodeNames[device.terminals[1]]
    n2 = nodeNames[device.terminals[2]]
    A[n1, n1] *= "+1/"*device.name*".R"
    A[n1, n2] *= "-1/"*device.name*".R"
    A[n2, n1] *= "-1/"*device.name*".R"
    A[n2, n2] *= "+1/"*device.name*".R"
end

function RHSSym(A::Vector, device::Resistor, nodeNames::Dict{String,Int}, lastSolution::Vector,  frequency::Float64=0.0)
    #Do nothing
end

function amountOfMatrixEntries(device::Resistor)
    return 2
end

function nameOfMatrixEntries(device::Resistor)
    return device.terminals
end