using Parameters
using LinearAlgebra
using Plots
using CSV
using DataFrames
include("devices/resistor.jl")
include("devices/voltagesource.jl")
include("devices/currentsource.jl")
include("devices/diode.jl")
include("devices/opv.jl")
include("devices/inductor.jl")
include("devices/capacitor.jl")
include("parsing.jl")


function assignDictNumber(dict::Dict{String, Int}, name)
    if !haskey(dict, name)
        dict[name] = length(dict)+1
    end
end

"Assigns an integer number to each node occuring in the circuit, starting with 1 for ground (which needs to have name 0)"
function assignNodeAndVoltageNumbers(elements)
    nodeDict = Dict{String, Int}(
        "0" => 1
    )
    for element in values(elements)
        for terminal in nameOfMatrixEntries(element)
            assignDictNumber(nodeDict, terminal)
        end
    end
    return nodeDict
end

"Lets all elements write their admittances to the matrix and rhs, skips one row and column for ground (otherwise we would be overdetermined)"
function assembleMatrixAndRhs(elements :: Dict, nodeDict::Dict{String, Int}, type::Type)
    A = zeros(type, length(nodeDict), length(nodeDict))
    rhs = zeros(type, length(nodeDict))
    for element in values(elements)
        writeToMatrix(A, element, nodeDict, [], 0.0)
        RHS(rhs, element, nodeDict, [], 0.0)
    end
    return A[2:end,2:end], rhs[2:end]
end

"Lets all elements write their admittances to the matrix and rhs, skips one row and column for ground (otherwise we would be overdetermined)"
function assembleMatrixAndRhs(elements :: Dict, nodeDict::Dict{String, Int}, type::Type, lastSolution=[])
    A = zeros(type, length(nodeDict), length(nodeDict))
    rhs = zeros(type, length(nodeDict))
    for element in values(elements)
        writeToMatrix(A, element, nodeDict, lastSolution, 0.0)
        RHS(rhs, element, nodeDict, lastSolution, 0.0)
    end
    return A[2:end,2:end], rhs[2:end]
end

"Lets all elements write their admittances to the matrix and rhs, skips one row and column for ground (otherwise we would be overdetermined)"
function assembleMatrixAndRhs(elements :: Dict, nodeDict::Dict{String, Int},type::Type, frequency::Real, lastSolution=[])
    A = zeros(type, length(nodeDict), length(nodeDict))
    rhs = zeros(type, length(nodeDict))
    for element in values(elements)
        writeToMatrix(A, element, nodeDict, lastSolution, frequency)
        RHS(rhs, element, nodeDict, lastSolution, frequency)
    end
    return A[2:end,2:end], rhs[2:end]
end


function solveEquation(elements, nodeDict,type::Type = Float64,frequency=0.0, lastSolution=zeros(type, length(nodeDict)))
    A,rhs = assembleMatrixAndRhs(elements, nodeDict, type)
    Ainv = inv(A)
    v = Ainv*rhs
    return [0;v]
end

function solveEquationNonLinear(elements, nodeDict, lastSolution=zeros(Real, length(nodeDict)))
    converged = false
    iterations = 0
    while converged == false
        A,rhs = assembleMatrixAndRhs(elements, nodeDict, lastSolution=lastSolution)
        v = A\rhs
        if(norm([0;v] - lastSolution) < 1e-12)
            converged = true
        end
        lastSolution = [0;v]
        iterations += 1
    end
    return lastSolution
end

"performs a DC sweep of voltage source V1 and returns node out"
function dcSweep(filename)
    elements = parseSpiceFile(filename)
    display(elements)
    nodeDict = assignNodeAndVoltageNumbers(elements)
    Vs = LinRange(5, 0, 10000)
    lastSolution = zeros(Real, length(nodeDict))
    res = Vector{Real}()
    for v in Vs
        elements["V1"].U = v
        lastSolution = solveEquation(elements, nodeDict)
        push!(res,lastSolution[nodeDict["out"]])
    end
    #data = CSV.read("output.csv", DataFrame, header=false, delim=' ', ignorerepeated=true)
    #plot(data[!,1], data[!,2])
    #plot!(Vector(Vs), res)
    #savefig("test.png")
    return Vs,res
end

"calculates the DC state including nonlinear components of some file"
function nonlinear(filename)
    elements = parseSpiceFile(filename)
    display(elements)
    nodeDict = assignNodeAndVoltageNumbers(elements)
    lastSolution = zeros(Real, length(nodeDict)-1)
    converged = false
    iterations = 0
    while converged == false
        A,rhs = assembleMatrixAndRhs(elements, nodeDict,Float64, [0;lastSolution])
        v = A\(rhs)
        println(iterations)
        println(norm(v - lastSolution))
        if(norm(v - lastSolution) < 1e-4)
            converged = true
        end
        lastSolution = v
        iterations += 1
    end
    display(lastSolution)
end




#This was a weird attemt to create a symbolic matrix, it works but is very slow
#Please mostly ignore this
function assembleMatrixAndRhsSym(elements :: Dict, nodeDict::Dict{String, Int})
    A = fill("0", length(nodeDict), length(nodeDict))
    rhs = fill("0", length(nodeDict))
    for element in values(elements)
        writeToMatrixSym(A, element, nodeDict, [], 0.0)
        RHSSym(rhs, element, nodeDict, [], 0.0)
    end
    return A[2:end,2:end], rhs[2:end]
end

function assembleEnvironment(elements :: Dict)
    env = String[]
    for element in values(elements)
        push!(env, element.name*"="*repr(element))
    end
    return join(env, ";")
end

function printResults(v, nodes)
    
    for n in keys(nodes)
        num = nodes[n]
        println("Node ", n, ": ", num , ' ' , v[num])
    end
    
end

function symbolicMatrix()
    elements = parseSpiceFile("test.cir")
    nodeDict = assignNodeAndVoltageNumbers(elements)
    mat,rhs = assembleMatrixAndRhsSym(elements, nodeDict)
    env = assembleEnvironment(elements)
    eval(Meta.parse(env))
    print(R1)
    matr = Meta.parse.(mat)
    Rhs = Meta.parse.(rhs)
    Vs = LinRange(5, 0, 1000000)
    res = []
    print(Rhs)
    for v in Vs
        V1.U = v
        s = eval.(matr)\eval.(Rhs)
        push!(res,s[nodeDict["2"]-1])
    end
    #plot(Vs, res)
end

