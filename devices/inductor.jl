mutable struct Inductor
    L
    terminals::Vector{String}
    name::String
end

function writeToMatrix(A::Matrix, device::Inductor, nodeNames::Dict{String,Int},  lastSolution::Vector, frequency::Float64=0.0) 
    n1 = nodeNames[device.terminals[1]]
    n2 = nodeNames[device.terminals[2]]
    admittance = 0.0
    if(frequency != 0)
        admittance = 1/(device.L * frequency * 1im)
    end
    A[n1, n1] += admittance
    A[n1, n2] -= admittance
    A[n2, n1] -= admittance
    A[n2, n2] += admittance
end

function RHS(A::Vector, device::Inductor, nodeNames::Dict{String,Int}, lastSolution::Vector, frequency::Float64=0.0)
    #Do nothing
end


function nameOfMatrixEntries(device::Inductor)
    return device.terminals
end