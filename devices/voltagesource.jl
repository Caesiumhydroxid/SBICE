mutable struct VoltageSource
    U
    terminals::Vector{String}
    name::String
end

function writeToMatrix(A::Matrix, device::VoltageSource, nodeNames::Dict{String,Int},  lastSolution::Vector, frequency::Float64=0.0) 
    voltageEntry = nodeNames[device.name] 
    n1 = nodeNames[device.terminals[1]]
    n2 = nodeNames[device.terminals[2]]
    A[n1, voltageEntry] += 1
    A[n2, voltageEntry] -= 1
    A[voltageEntry, n1] += 1
    A[voltageEntry, n2] -= 1
end

function RHS(A::Vector, device::VoltageSource, nodeNames::Dict{String,Int}, lastSolution::Vector,  frequency::Float64=0.0)
    voltageEntry = nodeNames[device.name] 
    A[voltageEntry] -= device.U
end

function writeToMatrixSym(A::Matrix, device::VoltageSource, nodeNames::Dict{String,Int},  lastSolution::Vector, frequency::Float64=0.0) 
    voltageEntry = nodeNames[device.name] 
    n1 = nodeNames[device.terminals[1]]
    n2 = nodeNames[device.terminals[2]]
    A[n1, voltageEntry] *= "+1"
    A[n2, voltageEntry] *= "-1"
    A[voltageEntry, n1] *= "+1"
    A[voltageEntry, n2] *= "-1"
end

function RHSSym(A::Vector, device::VoltageSource, nodeNames::Dict{String,Int}, lastSolution::Vector,  frequency::Float64=0.0)
    voltageEntry = nodeNames[device.name] 
    A[voltageEntry] *= "-"*device.name*".U"
end

function amountOfMatrixEntries(device::VoltageSource)
    return 3
end

function nameOfMatrixEntries(device::VoltageSource)
    return [device.terminals;device.name]
end