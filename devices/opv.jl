mutable struct OPV
    A
    terminals::Vector{String}
    name::String
end

function writeToMatrix(A::Matrix, device::OPV, nodeNames::Dict{String,Int},  lastSolution::Vector, frequency::Float64=0.0) 
    p = nodeNames[device.terminals[1]]
    n = nodeNames[device.terminals[2]]
    op = nodeNames[device.terminals[3]]
    on = nodeNames[device.terminals[4]]
    v = nodeNames[device.name]
    A[op,v] += 1
    A[on,v] -= 1
    A[v,op] += 1
    A[v,on] -= 1
    A[v,p] -= device.A
    A[v,n] += device.A
    
end

function RHS(A::Vector, device::OPV, nodeNames::Dict{String,Int}, lastSolution::Vector, frequency::Float64=0.0)
    #Do nothing
end

function amountOfMatrixEntries(device::OPV)
    return 5
end

function nameOfMatrixEntries(device::OPV)
    return [device.terminals;device.name]
end