using Turing
using StatsPlots
using ReverseDiff
pyplot()
include("../devices.jl")
include("../parsing.jl")
include("../probabilistic.jl")

tolerances = parseToleranceFile("measurements/spannungsteiler/rrr.tol")
measurementData = readMeasurementDataDC("measurements/spannungsteiler/RRR copy.CSV")
measurementData = measurementData[begin+10000:500:end-10000,:]
display(measurementData[!,"in"])
println(measurementData)
plot(measurementData[!,"in"],measurementData[!,"2"],label="ϕ₂")
plot!(measurementData[!,"in"],measurementData[!,"3"],label="ϕ₃")
xlabel!("Uin [V]")
ylabel!("Node Potentials[V]")
savefig("divider_data.pdf")
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

chn = sample(dcanalysis_uniform_prior(elements,nodes,voltages,amountOfMeasurements,length(observedNodeIndices),measurementMatrix,observedNodeIndices), NUTS(0.65), 10000)
test = chn[:,[Symbol("elementsR[\"R1\"].R"),Symbol("elementsR[\"R2\"].R"),Symbol("elementsR[\"R3\"].R")],:]
plot(test, background_color=:transparent)
savefig("divider_sampling.svg")

test = chn[:,[Symbol("offsets[1]"),Symbol("offsets[2]"),Symbol("noise[1]"),Symbol("noise[2]")],:]
plot(test,background_color=:transparent)
savefig("divider_sampling_more.svg")

R1 = vec(chn[Symbol("elementsR[\"R1\"].R")].data)
R2 = vec(chn[Symbol("elementsR[\"R2\"].R")].data)
hexbin(R1,R2,nbins=50, background_color=:transparent)
xlabel!("R1 [Ω]")
ylabel!("R2 [Ω]")
savefig("divider_sampling_hexbin.svg")

using Serialization
f = open("measurements/spannungsteiler/divider_chain_uniform1.jls", "w")
serialize(f,chn)
close(f)
weights = exp.(chn[:lp]);
p = weights / sum(weights);

