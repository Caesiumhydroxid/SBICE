using Turing
using StatsPlots
using ReverseDiff
include("devices.jl")
include("parsing.jl")
include("probabilistic.jl")

tolerances = parseToleranceFile("measurements/spannungsteiler/rrr.tol")
measurementData = readMeasurementDataDC("measurements/spannungsteiler/RRR copy.CSV")
measurementData = measurementData[begin+10000:500:end-10000,:]
display(measurementData[!,"in"])
println(measurementData)
plot(measurementData[!,"in"],measurementData[!,"2"])
plot!(measurementData[!,"in"],measurementData[!,"3"])

elements = parseSpiceFile("measurements/spannungsteiler/rrr.cir")
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

#missMeasurements = deepcopy(measurementData)
#allowmissing!(missMeasurements)
#missMeasurements[:,"2"] .= missing
#missMeasurements[:,"3"] .= missing
#display(missMeasurements)
#voltages = vec(measurementData[!,"in"])
#measurementMatrix=Matrix(select(measurementData, Not([:in, :t])))
#amountOfMeasurements = size(measurementMatrix,1)
#model = gdemo(elements,nodes,voltages,amountOfMeasurements,length(observedNodeIndices),measurementMatrix,observedNodeIndices)
#f,res,els = model()

voltages = vec(measurementData[!,"in"])
measurementMatrix=Matrix(select(measurementData, Not([:in, :t])))
amountOfMeasurements = size(measurementMatrix,1)
chn = sample(gdemo(elements,nodes,voltages,amountOfMeasurements,length(observedNodeIndices),measurementMatrix,observedNodeIndices), NUTS(0.65), 5000)
test = chn[:,[Symbol("elementsR[\"R1\"].R"),Symbol("elementsR[\"R2\"].R"),Symbol("elementsR[\"R3\"].R")],:]
plot(chn)
savefig("gdemo.png")

R1 = vec(chn[Symbol("elementsR[\"R1\"].R")].data)
R2 = vec(chn[Symbol("elementsR[\"R2\"].R")].data)
hexbin(R1,R2,nbins=50)

weights = exp.(chn[:lp]);
p = weights / sum(weights);

using StatsBase

Rs = vec(chn[Symbol("elementsR[\"R2\"].R")].data)
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