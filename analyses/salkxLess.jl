using DynamicHMC
using Turing
using StatsPlots
using ReverseDiff
include("../devices.jl")
include("../parsing.jl")
include("../probabilistic.jl")

tolerances = parseToleranceFile("measurements/sallenkeyExt/sallenkeyx.tol")
measuredNodes = ["1","2","out"]
measurementData = readMeasurementDataAC(["measurements/sallenkeyExt/SALEX01.CSV","measurements/sallenkeyExt/SALEX02.CSV","measurements/sallenkeyExt/SALEXOUT.CSV"], measuredNodes, 1e4)
measurementData = measurementData[begin:5:end,:]


elements = parseSpiceFile("measurements/sallenkeyExt/sallenkeyx.cir")
nodes = assignNodeAndVoltageNumbers(elements)
observedNodeIndices = Vector{Integer}()
display(nodes)
for x in measuredNodes
    if x in keys(nodes)
        push!(observedNodeIndices,nodes[x]-1)
    else 
        error("Node $x not found in spice file")
    end
end
print(observedNodeIndices)




measurementErrors = [0.0001,0.01,0.01]

frequencies = vec(measurementData[!,"f"])
measurementMatrix=Matrix(select(measurementData, Not([:f])))
phaseMatrix = rad2deg.(angle.(measurementMatrix))
measurementMatrix = 20 .*log10.(abs.(measurementMatrix))

#display(elements)
amountOfMeasurements = size(measurementMatrix,1)
amountOfMeasurementLocations = length(observedNodeIndices)

#elements["R1"].R = 100e3
#elements["R2"].R = 99.74e3
#elements["R3"].R = 149e3
#elements["R4"].R = 99.95e3
#elements["C1"].C = 2.70e-9
#elements["C2"].C = 3.22e-9

#tolerances["R1"] = 0.01
#tolerances["R2"] = 0.01
#tolerances["R3"] = 0.01
#tolerances["R4"] = 0.01
#tolerances["C1"] = 0.01
#tolerances["C2"] = 0.01

#model = acanalysis(elements,nodes,frequencies,amountOfMeasurements, amountOfMeasurementLocations, missing , missing,observedNodeIndices, measurementErrors)
#f,res,mes,phs,els = model()
#plot(f,phs[:,2],xscale=:log10)
#plot!(f,phaseMatrix[:,5],xscale=:log10)
#plot(f,mes[:,3],xscale=:log10)
#plot!(f,measurementMatrix[:,3],xscale=:log10)
#plot!(measurementData[!,"f"],20 .*log10.(abs.(measurementData[!,"3"])),xscale=:log10)

chn = sample(acanalysis(elements,nodes,frequencies,amountOfMeasurements, amountOfMeasurementLocations, measurementMatrix,phaseMatrix, observedNodeIndices, measurementErrors), Turing.NUTS(), 10_000)
test = chn[:,[Symbol("elementsR[\"R1\"].R"),Symbol("elementsR[\"C1\"].C"),Symbol("elementsR[\"R3\"].R")],:]
plot(chn)
using Serialization
f = open("measurements/sallenkeyExt/divider_chain_less.jls", "w")
serialize(f,chn)
close(f)

using PairPlots
using CairoMakie
chn2 = sample(acanalysis(elements,nodes,frequencies,amountOfMeasurements, amountOfMeasurementLocations, measurementMatrix,observedNodeIndices, measurementErrors), Turing.Prior(), 10_000)
fig = pairplot(chn,chn2)
save("pairplot_sallenkey_vs_prior.svg",fig)
fig = pairplot(chn)
save("pairplot_sallenkey.svg",fig)

plot(plots,layout=(length(variables),length(variables)),background_color=:transparent)
