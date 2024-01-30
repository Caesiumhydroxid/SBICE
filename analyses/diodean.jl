using Turing
using StatsPlots
using Plots
using ReverseDiff
pyplot()
include("../devices.jl")
include("../parsing.jl")
include("../probabilistic.jl")

## Parse circuit and measurement, cut out interesting measurement
tolerances = parseToleranceFile("measurements/diode_1n4148/diode.tol")
measurementData = readMeasurementDataDC("measurements/diode_1n4148/4148-01.CSV")
measurementData = measurementData[begin+12000:500:end-15000,:]
plot(measurementData[!,"in"],measurementData[!,"1"],label="ϕ₂")
plot!(measurementData[!,"in"],measurementData[!,"2"],label="ϕ₃")
xlabel!("Uin [V]")
ylabel!("Node Potentials[V]")
elements = parseSpiceFile("measurements/diode_1n4148/diode.cir")
nodes = assignNodeAndVoltageNumbers(elements)

observedNodeIndices = Vector{Integer}()
measuredNodes = intersect(setdiff(keys(nodes),["0","in"]),names(measurementData))
print(measuredNodes)
for x in keys(nodes)
    if x in measuredNodes
        push!(observedNodeIndices,nodes[x]-1)
    end
end

voltages = vec(measurementData[!,"in"])
measurementMatrix=Matrix(select(measurementData, Not([:in, :t])))
amountOfMeasurements = size(measurementMatrix,1)

## Perform Inference
chn = sample(dcanalysisnonlinear(elements,nodes,voltages,amountOfMeasurements,length(observedNodeIndices),measurementMatrix,observedNodeIndices), NUTS(0.65), 10000)

using Serialization
f = open("measurements/diode_1n4148/diode_1.jls", "w")
serialize(f,chn)
close(f)
chn = deserialize("measurements/diode_1n4148/diode_1.jls")


# Additional model for creating visualizations
@model function dcanalysisnonlinearSim(elements, nodes, voltages::Vector{Float64},amountOfMeasurements::Integer, observedNodeIndicesLength::Integer, measurementMatrix, observedNodeIndices::Vector{Integer}, noiseSim, offsetSim, ::Type{T} = Float64) where {T}
    elementsR = deepcopy(elements)
    type = T
    for e in keys(tolerances)
        if e[1] == 'R'
            elementsR[e].R = missing
            elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e]) 
            type = typeof(elementsR[e].R)
        end
        if e[1] == 'D'
            res = split(e,".")
            #check if string ends with IS
            if res[2] == "IS"
                elementsR[res[1]].IS = missing
                elementsR[res[1]].IS ~ Normal(elements[res[1]].IS*1e9, elements[res[1]].IS *1e9* tolerances[e])
                elementsR[res[1]].IS = elementsR[res[1]].IS * 1e-9
            end
            if res[2] == "N"
                elementsR[res[1]].N = missing
                elementsR[res[1]].N ~ Normal(elements[res[1]].N, elements[res[1]].N * tolerances[e])
            end
            if res[2] == "Gmin"
                elementsR[res[1]].Gmin = missing
                elementsR[res[1]].Gmin ~ Normal(elements[res[1]].Gmin, elements[res[1]].Gmin * tolerances[e])
            end
        end
    end

    offsets = offsetSim
    noise = noiseSim

    
    results = zeros(type,amountOfMeasurements,observedNodeIndicesLength)
    for i in 1:amountOfMeasurements
        lastSolution = zeros(type,length(nodes)-1)
        iterations = 0
        converged = false
        while converged == false
            A,rhs = assembleMatrixAndRhs(elementsR,nodes,type,[0;lastSolution])
            rhs[nodes[elementsR["V1"].name]-1] = -voltages[i]
            v = A\(rhs)
            if(norm(v - lastSolution) < 1e-5)
                converged = true
            end
            if(iterations > 1000)
                error("No convergence")
            end
            lastSolution = v
            iterations += 1
        end
        results[i,:] = lastSolution[observedNodeIndices] 
    end
    
    if measurementMatrix === missing
        measurementMatrix = zeros(type,amountOfMeasurements,observedNodeIndicesLength)
    end
    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(col.+offsets[i], diagm(noise[i]*ones(amountOfMeasurements)))
    end
    return voltages,measurementMatrix
end

tolerances = parseToleranceFile("measurements/diode_1n4148/diode.tol")
elements = parseSpiceFile("measurements/diode_1n4148/diode.cir")
model = dcanalysisnonlinearSim(elements,nodes,voltages,amountOfMeasurements,length(observedNodeIndices),missing,observedNodeIndices,[1e-3,1e-3],[-0.002,-0.0035])
plot()
xlims!(0,2)
ylims!(0,1)
for i in 1:100
    voltages,measurements = model()
    plot!(voltages,measurements[:,1],label="Prior ϕ₂", linealpha=0.05)
    plot!(voltages,measurements[:,2],label="Prior ϕ₃", linealpha=0.05)
end
plot!(measurementData[!,"in"],measurementData[!,"1"],label="Measured ϕ₂")
plot!(measurementData[!,"in"],measurementData[!,"2"],label="Measured ϕ₃")
plot!(legend=false)
elements["D1"].IS = 16.68e-9
elements["R2"].R = 219
model = dcanalysisnonlinearSim(elements,nodes,voltages,amountOfMeasurements,length(observedNodeIndices),missing,observedNodeIndices,[1e-5,1e-5],[-0.002,-0.0035])
voltages,measurements = model()
plot!(voltages,measurements[:,1],label="Posterior ϕ₂", linestyle=:dot, linewidth = 2)
plot!(voltages,measurements[:,2],label="Posterior ϕ₃",linestyle=:dot, linewidth = 2)
plot!(background_color = :transparent)
savefig("diode_inference.svg")