using Turing
using StatsPlots
using ReverseDiff
include("../devices.jl")
include("../parsing.jl")
include("../probabilistic.jl")

tolerances = parseToleranceFile("measurements/sallenkey/sallenkey.tol")
measuredNodes = ["1","2","3","out"]
measurementData = readMeasurementDataAC(["measurements/sallenkey/SALEN1.CSV","measurements/sallenkey/SALEN2.CSV","measurements/sallenkey/SALEN3.CSV","measurements/sallenkey/SALENOUT.CSV"], measuredNodes, 1e4)
measurementData = measurementData[begin:5:end,:]
plot(measurementData[!,"f"] ,abs.(measurementData[!,"out"]),xscale=:log10,yscale=:log10)
plot(measurementData[!,"f"],abs.(measurementData[!,"3"]),xscale=:log10,yscale=:log10)

elements = parseSpiceFile("measurements/sallenkey/sallenkey.cir")
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




measurementErrors = [0.0001,0.001,0.001,0.001]

frequencies = vec(measurementData[!,"f"])
measurementMatrix=Matrix(select(measurementData, Not([:f])))
phaseMatrix = rad2deg.(angle.(measurementMatrix))
measurementMatrix = 20 .*log10.(abs.(measurementMatrix))

#display(elements)
amountOfMeasurements = size(measurementMatrix,1)
amountOfMeasurementLocations = length(observedNodeIndices)

model = acanalysis(elements,nodes,frequencies,amountOfMeasurements, amountOfMeasurementLocations, missing , missing,observedNodeIndices, measurementErrors)
f,res,mes,phs,els = model()
plot(f,phs[:,4],xscale=:log10)
plot(f,phaseMatrix[:,4],xscale=:log10)
#plot!(measurementData[!,"f"],20 .*log10.(abs.(measurementData[!,"3"])),xscale=:log10)
display(els)


chn = sample(acanalysis(elements,nodes,frequencies,amountOfMeasurements, amountOfMeasurementLocations, measurementMatrix, phaseMatrix,observedNodeIndices, measurementErrors), Turing.NUTS(), 10_000)
test = chn[:,[Symbol("elementsR[\"R1\"].R"),Symbol("elementsR[\"C1\"].C"),Symbol("elementsR[\"R3\"].R")],:]
plot(chn)
using Serialization
f = open("measurements/sallenkey/chain.jls", "w")
serialize(f,chn)
close(f)
savefig("gdemo.png")
print(tolerances)
R1 = vec(chn[Symbol("elementsR[\"R1\"].R")].data)
R2 = vec(chn[Symbol("elementsR[\"C1\"].C")].data)
hexbin(R1,R2,nbins=30)

weights = exp.(chn[:lp]);
p = weights / sum(weights);

using StatsBase

Rs = vec(chn[Symbol("elementsR[\"R1\"].R")].data)
h = fit(Histogram, Rs, nbins=50)
plot(h)
max_bin_index = argmax(h.weights)
bin_edges = h.edges[1][max_bin_index:max_bin_index+1]
estimated_max = sum(bin_edges) / 2

# Get maximum a posteriori (MAP) estimate
histogram(Rs)
println(keys(chn))
map_estimate = chn[Symbol("elementsR[\"R1\"].R")][argmax(p)]
println(map_estimate)


map_estimate = mean(chn[Symbol("elementsR[\"R1\"].R")])
println(map_estimate)
map_estimate = mean(chn[Symbol("elementsR[\"R2\"].R")])
println(map_estimate)