using DynamicHMC
using Turing
using StatsPlots
using ReverseDiff
include("../devices.jl")
include("../parsing.jl")
include("../probabilistic.jl")

tolerances = parseToleranceFile("measurements/sallenkeyExt/sallenkeyx.tol")
measuredNodes = ["1","2","3","4","out"]
measurementData = readMeasurementDataAC(["measurements/sallenkeyExt/SALEX01.CSV","measurements/sallenkeyExt/SALEX02.CSV","measurements/sallenkeyExt/SALEX03.CSV","measurements/sallenkeyExt/SALEX04.CSV","measurements/sallenkeyExt/SALEXOUT.CSV"], measuredNodes, 1e4)
measurementData = measurementData[begin:5:end,:]
A1 = plot(measurementData[!,"f"] ,abs.(measurementData[!,"out"]),xscale=:log10,yscale=:log10,label="ϕ₅")
A1 = plot!(measurementData[!,"f"],abs.(measurementData[!,"3"]),xscale=:log10,yscale=:log10,label="ϕ₃")
A1 = plot!(measurementData[!,"f"],abs.(measurementData[!,"2"]),xscale=:log10,yscale=:log10,label="ϕ₂")
A1 = plot!(measurementData[!,"f"],abs.(measurementData[!,"1"]),xscale=:log10,yscale=:log10,label="ϕ₁")
A1 = plot!(measurementData[!,"f"],abs.(measurementData[!,"4"]),xscale=:log10,yscale=:log10,label="ϕ₄")
A1 = xaxis!("f [Hz]")
A1 = yaxis!("|U| [V]")
A2 = plot(measurementData[!,"f"] ,(180/pi).*(angle.(measurementData[!,"out"])),xscale=:log10,label="ϕ₅")
A2 = plot!(measurementData[!,"f"] ,(180/pi).*(angle.(measurementData[!,"1"])),xscale=:log10,label="ϕ₁")
A2 = plot!(measurementData[!,"f"] ,(180/pi).*(angle.(measurementData[!,"2"])),xscale=:log10,label="ϕ₂")
A2 = plot!(measurementData[!,"f"] ,(180/pi).*(angle.(measurementData[!,"3"])),xscale=:log10,label="ϕ₃")
A2 = plot!(measurementData[!,"f"] ,(180/pi).*(angle.(measurementData[!,"4"])),xscale=:log10,label="ϕ₄")
A2 = xaxis!("f [Hz]")
A2 = yaxis!("ϕ [°]")
plot(A1,A2,layout=(2,1), background_color=:transparent)
savefig("sallenkeyx_data.svg")

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




measurementErrors = [0.0001,0.001,0.001,0.001,0.001]

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
#plot(f,phs[:,5],xscale=:log10)
#plot!(f,phaseMatrix[:,5],xscale=:log10)
#plot(f,mes[:,4],xscale=:log10)
#plot!(f,measurementMatrix[:,4],xscale=:log10)
#plot!(measurementData[!,"f"],20 .*log10.(abs.(measurementData[!,"3"])),xscale=:log10)

chn = sample(acanalysis(elements,nodes,frequencies,amountOfMeasurements, amountOfMeasurementLocations, measurementMatrix,observedNodeIndices, measurementErrors), Turing.NUTS(), 10_000)
test = chn[:,[Symbol("elementsR[\"R1\"].R"),Symbol("elementsR[\"C1\"].C"),Symbol("elementsR[\"R3\"].R")],:]
plot(chn)
using Serialization
f = open("measurements/sallenkeyExt/divider_chain.jls", "w")
serialize(f,chn)
close(f)

using PairPlots
using CairoMakie
chn2 = sample(acanalysis(elements,nodes,frequencies,amountOfMeasurements, amountOfMeasurementLocations, measurementMatrix,observedNodeIndices, measurementErrors), Turing.Prior(), 10_000)
pairplot(chn,chn2)
savefig("pairplot_sallenkey_vs_prior.svg")

variables = chn.name_map.parameters[[1,2,3,4,5,7]]
plots = Matrix{Plots.Plot}(undef, length(variables),length(variables))
for (i,v) in enumerate(variables)
    for (ii,vv) in enumerate(variables)
        if v != vv
            v1 = vec(chn[v].data)
            v2 = vec(chn[vv].data)
            plots[i,ii] = hexbin(v1, v2, background_color=:transparent)
        end
    end
end
plot(plots,layout=(length(variables),length(variables)),background_color=:transparent)

savefig("gdemo.png")
print(tolerances)
R1 = vec(chn[Symbol("elementsR[\"R3\"].R")].data)
R2 = vec(chn[Symbol("elementsR[\"R4\"].R")].data)
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