using CSV

function parseSIValue(string)
    print(string)
    res = match(r"([\d.e]+)([a-zA-Z]?)", string)
    print(res)
    value = parse(Float64, res[1])
    
    if(res[2] == "p")
        value *= 1e-12
    elseif(res[2] == "n")
        value *= 1e-9
    elseif(res[2] == "u")
        value *= 1e-6
    elseif(res[2] == "m")
        value *= 1e-3
    elseif(res[2] == "k")
        value *= 1e3
    elseif(res[2] == "X")
        value *= 1e6
    elseif(res[2] == "G")
        value *= 1e9
    elseif(res[2] == "T")
        value *= 1e12
    end
    return value
end

function parseModel(models::Dict{String,Any}, line)
    #Match regex
    result = match(r".MODEL\s(\w+)\s(\w)\s+\(([^\)]+)\)", line)
    if !isnothing(result)
        modelName = result.captures[1]
        modelType = result.captures[2] 
        modelParams = result.captures[3]
        modelParams = split(modelParams, " ")
        filter!(x -> x != "", modelParams)
        println(modelParams)
        if(modelType == "D")
            diode = Diode(name=modelName)
            for param in modelParams
                pair = split(param, "=")
                key = strip(pair[1])
                value = parseSIValue(strip(pair[2]))
                println(key,typeof(value),value)
                setfield!(diode, Symbol(key), value)
            end
            models[modelName] = diode
        end
    end
end

function parseSpiceFile(inputFile::String)

    #Read the file line by line
    lines = readlines(inputFile)
    newLines = Vector{String}()
    #Merge lines that start with a '+' with the previous line
    for i in eachindex(lines)
        line = lines[i]*' '
        if line[1] == '+'
            newLines[length(newLines)] = newLines[length(newLines)] * line[2:end]
        else
            push!(newLines, line)
        end
    end

    elements = Dict{String,Any}()
    models = Dict{String, Any}()
    ##Parse elements into list
    for line in newLines
        println(line)
        if line[1] == '.'
            parseModel(models,line)
        end
        line=strip(line)
        if(length(line)==0)
            continue
        end
        if line[1] == 'R'
            result = split(line, " ")
            name = result[1]
            node1 = result[2]
            node2 = result[3]
            resistance = parseSIValue(result[4])
            elements[name]=Resistor(resistance, [node1, node2], name)
        elseif line[1] == 'C'
            result = split(line, " ")
            name = result[1]
            node1 = result[2]
            node2 = result[3]
            capacitance = parseSIValue(result[4])
            elements[name]=Capacitor(capacitance, [node1, node2], name)
        elseif line[1] == 'L'
            result = split(line, " ")
            name = result[1]
            node1 = result[2]
            node2 = result[3]
            inductance = parseSIValue(result[4])
            elements[name]=Inductor(inductance, [node1, node2], name)
        elseif line[1] == 'I'
            result = split(line, " ")
            name = result[1]
            node1 = result[2]
            node2 = result[3]
            type = result[4]
            current = parseSIValue(result[4])
            elements[name]=CurrentSource(current, [node1, node2], name)
        elseif line[1] == 'V'
            result = split(line, " ")
            name = result[1]
            node1 = result[2]
            node2 = result[3]
            type = result[4]
            voltage = parseSIValue(result[5])
            elements[name]=VoltageSource(voltage, [node1, node2], name)
        elseif line[1] == 'D'
            result = split(line, " ")
            name = result[1]
            node1 = result[2]
            node2 = result[3]
            model = result[4]
            diode = copy(models[model])
            diode.name = name
            diode.terminals = [node1, node2]
            init(diode)
            elements[name]=diode
        elseif line[1]== 'O'
            result = split(line, " ")
            name = result[1]
            np = result[2]
            nn = result[3]
            npo = result[4]
            nno = result[5]
            gain = parseSIValue(result[6])
            terminals = [np,nn,npo,nno]
            elements[name]=OPV(gain, terminals, name)
        else
            println("Unknown element: ", line)
        end

    end
    return elements
end

function parseToleranceNumber(text)
    text = strip(text)
    if text[end] == '%'
        return parse(Float64, text[1:end-1])/100    
    else
        return parse(Float64, text)
    end
end

function parseToleranceFile(inputFile::String)
    lines = readlines(inputFile)
    tolerances = Dict{String,Any}()
    for line in lines
        result = split(line, " ")
        name = result[1]
        tolerance = parseToleranceNumber(result[2])
        tolerances[name]=tolerance
    end
    return tolerances
end


function readMeasurementDataDC(inputFile::String)
    data = CSV.read(inputFile,DataFrame, delim=',', header=1)
    return data
end

function readMeasurementDataAC(inputFiles, nodes, fmax=1e5)
    nodenames = ["f"; nodes]
    data = nothing
    first = true
    for (file,node) in zip(inputFiles,nodes)
        csv = CSV.read(file,DataFrame, delim=',', header=1)
        if(first)
            data = DataFrame([csv[!, Symbol("Frequency in Hz")] for _ = nodenames] , nodenames)
            first = false
        end
        display(csv)
        data[!, node] = (10 .^(csv[!, Symbol("Gain in dB")]/20)) .* (cos.(csv[!, Symbol("Phase in °")]*pi/180) + im*sin.(csv[!, Symbol("Phase in °")]*pi/180))
    end
    return data[(data.f .<= fmax), :]
end


function parseTolerances(inputFile::String)
    lines = readlines(inputFile)
    elements = Dict{String,Any}()
    for line in lines
        if line[1] == 'R'
            result = split(line, " ")
            name = result[1]
            tolerance = parse(Float64,result[2])
            elements[name]=tolerance
        end
        if line[1] == 'C'
            result = split(line, " ")
            name = result[1]
            tolerance = parse(Float64,result[2])
            elements[name]=tolerance
        end
        if line[1] == 'L'
            result = split(line, " ")
            name = result[1]
            tolerance = parse(Float64,result[2])
            elements[name]=tolerance
        end
        if line[1] == 'I'
            result = split(line, " ")
            name = result[1]
            tolerance = parse(Float64,result[2])
            elements[name]=tolerance
        end
        if line[1] == 'V'
            result = split(line, " ")
            name = result[1]
            tolerance = parse(Float64,result[2])
            elements[name]=tolerance
        end
        if line[1] == 'D'
            result = split(line, " ")
            name = result[1]
            tolerance = parse(Float64,result[2])
            elements[name]=tolerance
        end
    end
end