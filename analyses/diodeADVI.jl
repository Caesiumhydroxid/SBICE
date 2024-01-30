using Turing
using StatsPlots
using Plots
using ReverseDiff
pyplot()
include("../devices.jl")
include("../parsing.jl")
include("../probabilistic.jl")

tolerances = parseToleranceFile("measurements/diode_1n4148/diode.tol")
#measurementData = readMeasurementDataDC("measurements/diode_1n4148/4148-04.CSV")
measurementData = readMeasurementDataDC("measurements/diode_4004/4004.CSV")
measurementData = measurementData[begin+12000:500:end-15000,:]
plot(measurementData[!,"in"],measurementData[!,"1"],label="ϕ₂")
plot!(measurementData[!,"in"],measurementData[!,"2"],label="ϕ₃")
xlabel!("Uin [V]")
ylabel!("Node Potentials[V]")
elements = parseSpiceFile("measurements/diode_1n4148/diode.cir")
nodes = assignNodeAndVoltageNumbers(elements)

measurements = measurementData
observedNodeIndices = Vector{Integer}()
measuredNodes = intersect(setdiff(keys(nodes),["0","in"]),names(measurements))
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
model = dcanalysisnonlinear(elements,nodes,voltages,amountOfMeasurements,length(observedNodeIndices),measurementMatrix,observedNodeIndices)
using Optim
#This tends to fail because we get a NAN or Inf in the matrix. Just try a bunch of times :)
map_estimate = optimize(model,MAP())
display(map_estimate.values)

using StatsBase
coeftable(map_estimate)