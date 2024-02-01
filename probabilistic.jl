using Turing
include("devices.jl")

@model function dcanalysis(elements, nodes, voltages::Vector{Float64},amountOfMeasurements::Integer, observedNodeIndicesLength::Integer, measurementMatrix::Matrix{Float64}, observedNodeIndices::Vector{Integer}, ::Type{T} = Float64) where {T}   
    elementsR = deepcopy(elements)
    type = T
    for e in keys(tolerances)
        if e[1] == 'R'
            elementsR[e].R = missing
            elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e]) 
            type = typeof(elementsR[e].R)
        end
    end

    offsets ~ MvNormal(zeros(observedNodeIndicesLength), diagm(0.01.*ones(observedNodeIndicesLength)))
    noise = zeros(observedNodeIndicesLength)
    for i in 1:observedNodeIndicesLength
        noise[i]   ~ Uniform(0, 0.01)
    end
    
    A,RHS = assembleMatrixAndRhs(elementsR,nodes,type)
    Ainv = inv(A)
    
    results = zeros(type,amountOfMeasurements,observedNodeIndicesLength)
    for i in 1:amountOfMeasurements
        RHS[nodes[elementsR["V1"].name]-1] = -voltages[i]
        res = Ainv*RHS
        results[i,:] = res[observedNodeIndices] 
    end

    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(col.+offsets[i], diagm(noise[i]*ones(amountOfMeasurements)))
    end
    return voltages,measurementMatrix,elementsR
end

@model function dcanalysis_uniform_prior(elements, nodes, voltages::Vector{Float64},amountOfMeasurements::Integer, observedNodeIndicesLength::Integer, measurementMatrix::Matrix{Float64}, observedNodeIndices::Vector{Integer}, ::Type{T} = Float64) where {T}
    
    elementsR = deepcopy(elements)
    type = T
    for e in keys(tolerances)
        if e[1] == 'R'
            elementsR[e].R = missing
            elementsR[e].R ~ Uniform(elements[e].R * (1-tolerances[e]), elements[e].R * (1+tolerances[e])) 
            type = typeof(elementsR[e].R)
        end
    end

    noise = zeros(observedNodeIndicesLength)
    offsets = zeros(observedNodeIndicesLength)
    for i in 1:observedNodeIndicesLength
        noise[i]   ~ Uniform(0, 0.01)
        offsets[i] ~ Uniform(-0.05, 0.05)
    end
    
    A,RHS = assembleMatrixAndRhs(elementsR,nodes,type)
    Ainv = inv(A)
    
    results = zeros(type,amountOfMeasurements,observedNodeIndicesLength)
    for i in 1:amountOfMeasurements
        RHS[nodes[elementsR["V1"].name]-1] = -voltages[i]
        res = Ainv*RHS
        results[i,:] = res[observedNodeIndices] 
    end

    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(col.+offsets[i], diagm(noise[i]*ones(amountOfMeasurements)))
    end
    return voltages,measurementMatrix,elementsR
end

@model function acanalysis(elements, nodes, frequencies::Vector{Float64},amountOfMeasurements::Integer, amountOfMeasurementLocations::Integer, measurementMatrix, observedNodeIndices::Vector{Integer}, measurementErrors, ::Type{T} = ComplexF64) where {T}
    elementsR = deepcopy(elements)
    type = Nothing
    for e in keys(tolerances)
        if e[1] == 'R'
            elementsR[e].R = missing
            #elementsR[e].R ~ truncated(Normal(elements[e].R, elements[e].R * tolerances[e]), elements[e].R * (1-5*tolerances[e]), elements[e].R * (1+5*tolerances[e]))
            elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e])
        end
        if e[1] == 'C'
            elementsR[e].C = missing
            elementsR[e].C ~ Normal(elements[e].C*1e9, elements[e].C *1e9 * tolerances[e]) 
            elementsR[e].C *= ComplexF64(1e-9)
            type = typeof(elementsR[e].C)
        end
        if e[1] == 'L'
            elementsR[e].L = missing
            elementsR[e].L ~ Normal(elements[e].L*1e6, elements[e].L*1e6 * tolerances[e]) 
            elementsR[e].L *= ComplexF64(1e-6)
            type = typeof(elementsR[e].L)
        end
    end
   
    #measurementErrors ~ filldist(InverseGaussian(1,0.001), amountOfMeasurementLocations)

    results = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    for (x,f) in enumerate(frequencies)
        A,RHS = assembleMatrixAndRhs(elementsR,nodes,type,f)
        res = A\RHS
        results[x,:] = res[observedNodeIndices] 
    end

    if measurementMatrix === missing
        measurementMatrix = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    end
    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(20 .*log10.(abs.(col)), diagm(measurementErrors[i]*ones(length(col))))
    end
    return frequencies,20 .*log10.(abs.(results)),measurementMatrix,elementsR
end

@model function acanalysis(elements, nodes, frequencies::Vector{Float64},amountOfMeasurements::Integer, amountOfMeasurementLocations::Integer, measurementMatrix, phaseMatrix, observedNodeIndices::Vector{Integer}, measurementErrors, ::Type{T} = ComplexF64) where {T}
    elementsR = deepcopy(elements)
    type = Nothing
    for e in keys(tolerances)
        if e[1] == 'R'
            elementsR[e].R = missing
            elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e]) 
        end
        if e[1] == 'C'
            elementsR[e].C = missing
            elementsR[e].C ~ Normal(elements[e].C*1e9, elements[e].C *1e9 * tolerances[e]) 
            elementsR[e].C *= ComplexF64(1e-9)
            type = typeof(elementsR[e].C)
        end
        if e[1] == 'L'
            elementsR[e].L = missing
            elementsR[e].L ~ Normal(elements[e].L*1e6, elements[e].L*1e6 * tolerances[e]) 
            elementsR[e].L *= ComplexF64(1e-6)
            type = typeof(elementsR[e].L)
        end
    end
   
    #measurementErrors ~ filldist(InverseGaussian(1,0.001), amountOfMeasurementLocations)

    results = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    for (x,f) in enumerate(frequencies)
        A,RHS = assembleMatrixAndRhs(elementsR,nodes,type,f)
        res = A\RHS
        results[x,:] = res[observedNodeIndices] 
    end
    if measurementMatrix === missing
        measurementMatrix = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    end
    if phaseMatrix === missing
        phaseMatrix = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    end
    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(20 .*log10.(abs.(col)), diagm(measurementErrors[i]*ones(length(col))))
        phaseMatrix[:,i] ~ MvNormal((180/pi).*(angle.(col)), diagm(0.1*ones(length(col))))
    end
    return frequencies,20 .*log10.(abs.(results)),measurementMatrix,phaseMatrix,elementsR
end

@model function acanalysisForVar(elements, nodes, frequencies::Vector{Float64},amountOfMeasurements::Integer, amountOfMeasurementLocations::Integer, measurementMatrix, phaseMatrix, observedNodeIndices::Vector{Integer}, ::Type{T} = ComplexF64) where {T}
    elementsR = deepcopy(elements)
    type = Nothing
    for e in keys(tolerances)
        if e[1] == 'R'
            elementsR[e].R = missing
            elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e]) 
        end
        if e[1] == 'C'
            elementsR[e].C = missing
            elementsR[e].C ~ Normal(elements[e].C*1e9, elements[e].C *1e9 * tolerances[e]) 
            elementsR[e].C *= ComplexF64(1e-9)
            type = typeof(elementsR[e].C)
        end
        if e[1] == 'L'
            elementsR[e].L = missing
            elementsR[e].L ~ Normal(elements[e].L*1e6, elements[e].L*1e6 * tolerances[e]) 
            elementsR[e].L *= ComplexF64(1e-6)
            type = typeof(elementsR[e].L)
        end
    end
   
    measurementErrors ~ filldist(InverseGaussian(1,0.001), amountOfMeasurementLocations)

    results = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    for (x,f) in enumerate(frequencies)
        A,RHS = assembleMatrixAndRhs(elementsR,nodes,type,f)
        res = A\RHS
        results[x,:] = res[observedNodeIndices] 
    end
    if measurementMatrix === missing
        measurementMatrix = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    end
    if phaseMatrix === missing
        phaseMatrix = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    end
    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(20 .*log10.(abs.(col)), diagm(measurementErrors[i]*ones(length(col))))
        phaseMatrix[:,i] ~ MvNormal((180/pi).*(angle.(col)), diagm(0.1*ones(length(col))))
    end
    return frequencies,20 .*log10.(abs.(results)),measurementMatrix,phaseMatrix,elementsR
end

@model function acanalysiswithfaults(elements, nodes, frequencies::Vector{Float64},amountOfMeasurements::Integer, amountOfMeasurementLocations::Integer, measurementMatrix, observedNodeIndices::Vector{Integer}, measurementErrors, ::Type{T} = ComplexF64) where {T}
    elementsR = deepcopy(elements)
    type = Nothing
    faultyElement = Dict{String,Int}()
    #shortOrOpen = Dict{String,Int}()
    for e in keys(tolerances)
        faultyElement[e] ~ Bernoulli(0.01)
        if e[1] == 'R'
            elementsR[e].R = missing
            if(faultyElement[e] == 1)
                elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e]) 
                elementsR[e] = Resistor(1e-8,elements[e].terminals,elements[e].name)
            else
               elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e]) 
               elementsR[e].R *= ComplexF64(1)
               type = typeof(elementsR[e].R)
            end
        end
        if e[1] == 'C'
            elementsR[e].C = missing
            if(faultyElement[e]==1)
                elementsR[e].C ~ Normal(elements[e].C*1e9, elements[e].C *1e9 * tolerances[e]) 
                elementsR[e] = Resistor(1e-8,elements[e].terminals,elements[e].name)
            else
                elementsR[e].C ~ Normal(elements[e].C*1e9, elements[e].C *1e9 * tolerances[e]) 
                elementsR[e].C *= ComplexF64(1e-9)
                type = typeof(elementsR[e].C)
            end
            
        end
    end

    results = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    for (x,f) in enumerate(frequencies)
        A,RHS = assembleMatrixAndRhs(elementsR,nodes,type,f)
        res = A\RHS
        results[x,:] = res[observedNodeIndices] 
    end

    if measurementMatrix === missing
        measurementMatrix = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    end
    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(20 .*log10.(abs.(col)), diagm(measurementErrors[i]*ones(length(col))))
    end
    return frequencies,20 .*log10.(abs.(results)),measurementMatrix,elementsR
end


@model function acanalysiswithfaultopen(elements, nodes, frequencies::Vector{Float64},amountOfMeasurements::Integer, amountOfMeasurementLocations::Integer, measurementMatrix, observedNodeIndices::Vector{Integer}, measurementErrors, ::Type{T} = ComplexF64) where {T}
    elementsR = deepcopy(elements)
    type = Nothing
    faultyElement = Dict{String,Int}()
    #shortOrOpen = Dict{String,Int}()
    for e in keys(tolerances)
        faultyElement[e] ~ Bernoulli(0.01)
        if e[1] == 'R'
            elementsR[e].R = missing
            if(faultyElement[e] == 1)
                elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e]) 
                elementsR[e] = Resistor(1e8,elements[e].terminals,elements[e].name)
            else
               elementsR[e].R ~ Normal(elements[e].R, elements[e].R * tolerances[e]) 
               elementsR[e].R *= ComplexF64(1)
               type = typeof(elementsR[e].R)
            end
        end
        if e[1] == 'C'
            elementsR[e].C = missing
            if(faultyElement[e]==1)
                elementsR[e].C ~ Normal(elements[e].C*1e9, elements[e].C *1e9 * tolerances[e]) 
                elementsR[e] = Resistor(1e8,elements[e].terminals,elements[e].name)
            else
                elementsR[e].C ~ Normal(elements[e].C*1e9, elements[e].C *1e9 * tolerances[e]) 
                elementsR[e].C *= ComplexF64(1e-9)
                type = typeof(elementsR[e].C)
            end
            
        end
    end

    results = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    for (x,f) in enumerate(frequencies)
        A,RHS = assembleMatrixAndRhs(elementsR,nodes,type,f)
        res = A\RHS
        results[x,:] = res[observedNodeIndices] 
    end

    if measurementMatrix === missing
        measurementMatrix = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    end
    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(20 .*log10.(abs.(col)), diagm(measurementErrors[i]*ones(length(col))))
    end
    return frequencies,20 .*log10.(abs.(results)),measurementMatrix,elementsR
end


@model function dcanalysisnonlinear(elements, nodes, voltages::Vector{Float64},amountOfMeasurements::Integer, observedNodeIndicesLength::Integer, measurementMatrix::Matrix{Float64}, observedNodeIndices::Vector{Integer}, ::Type{T} = Float64) where {T}
    
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

    offsets ~ MvNormal(zeros(observedNodeIndicesLength), diagm(0.01.*ones(observedNodeIndicesLength)))
    noise = zeros(observedNodeIndicesLength)
    for i in 1:observedNodeIndicesLength
        noise[i]   ~ Uniform(0, 0.01)
    end
    
    
    results = zeros(type,amountOfMeasurements,observedNodeIndicesLength)
    for i in 1:amountOfMeasurements
        lastSolution = zeros(type,length(nodes)-1)
        iterations = 0
        converged = false
        while converged == false
            A,rhs = assembleMatrixAndRhs(elementsR,nodes,type,[0;lastSolution])
            rhs[nodes[elementsR["V1"].name]-1] = -voltages[i]
            v = A\(rhs)
            if(norm(v - lastSolution) < 1e-4)
                converged = true
            end
            if(iterations > 10000)
                converged = true
            end
            lastSolution = v
            iterations += 1
        end
        results[i,:] = lastSolution[observedNodeIndices] 
    end
    
    if measurementMatrix === missing
        measurementMatrix = zeros(type,amountOfMeasurements,amountOfMeasurementLocations)
    end
    for (i,col) in enumerate(eachcol(results))
        measurementMatrix[:,i] ~ MvNormal(col.+offsets[i], diagm(noise[i]*ones(amountOfMeasurements)))
    end
    return voltages,measurementMatrix
end
