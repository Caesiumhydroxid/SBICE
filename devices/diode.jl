using Parameters
Base.@kwdef mutable struct Diode
    terminals :: Vector{String} = []
    name:: String = ""
    IS   = 10e-12
    RS   = 1e-12
    N    = 1
    TT   = 0
    CJO  = 0
    VJ   = 1
    M    = 0.5
    EG   = 1.11
    XTI  = 3
    KF   = 0
    AF   = 1
    FC   = 0.5
    BV   = Inf
    IBV  = 1e-3
    Gmin = 1e-12
end

Base.copy(s::Diode) = Diode(s.terminals, s.name, s.IS, s.RS, s.N, s.TT, s.CJO, s.VJ, s.M, s.EG, s.XTI, s.KF, s.AF, s.FC, s.BV, s.IBV, s.Gmin)


function init(device :: Diode)
    device.terminals = [device.terminals[1], device.terminals[1]*"primed",device.terminals[2]]

end

function writeToMatrix(A::Matrix, device::Diode, nodeNames::Dict{String,Int}, lastSolution::Vector, frequency::Float64=0.0) 
    n1 = nodeNames[device.terminals[1]]
    n1p = nodeNames[device.terminals[2]]
    n2 = nodeNames[device.terminals[3]]
    #Get the voltage across the diode
    vd = lastSolution[n1p] - lastSolution[n2]
    Vt = 0.026
    Vte = device.N * Vt
    
    if (vd >= -3 * Vte)
        evd = exp(vd / Vte);
        cd = device.IS * (evd - 1) + device.Gmin * vd;
        gd = device.IS * evd / Vte + device.Gmin;
    
    elseif (vd >= -device.BV)
        # Linear region
        arg = 3 * Vte / (vd * MathConstants.e);
        cd = -device.IS * (1 + arg^3) + device.Gmin * vd;
        gd = device.IS * 3 * (arg^3)/vd + device.Gmin;
    else
        evrev = exp(-(device.BV + vd) / Vte);
        cd = -device.IS * (evrev) + device.Gmin * vd;
        gd = device.IS * evrev / Vte + device.Gmin ;
    end

    #Write the current to the matrix
    A[n1, n1] += 1/device.RS
    A[n1, n1p] -= 1/device.RS
    A[n1p, n1] -= 1/device.RS
    A[n1p, n1p] += 1/device.RS
    A[n1p, n1p] += gd
    A[n1p, n2] -= gd
    A[n2, n1p] -= gd
    A[n2, n2] += gd
end

function RHS(A::Vector, device::Diode, nodeNames::Dict{String,Int}, lastSolution::Vector, frequency::Float64=0.0)
    n1 = nodeNames[device.terminals[1]]
    n1p = nodeNames[device.terminals[2]]
    n2 = nodeNames[device.terminals[3]]
    #Get the voltage across the diode
    vd = lastSolution[n1p] - lastSolution[n2]
    Vt = 0.026
    Vte = device.N * Vt
    
    if (vd >= -3 * Vte)
        evd = exp(vd / Vte);
        cd = device.IS * (evd - 1) + device.Gmin * vd;
        gd = device.IS * evd / Vte + device.Gmin;
    
    elseif (vd >= -device.BV)
        # Linear region
        arg = 3 * Vte / (vd * MathConstants.e);
        cd = -device.IS * (1 + arg^3) + device.Gmin * vd;
        gd = device.IS * 3 * (arg^3)/vd + device.Gmin;
    else
        evrev = exp(-(device.BV + vd) / Vte);
        cd = -device.IS * (evrev) + device.Gmin * vd;
        gd = device.IS * evrev / Vte + device.Gmin ;
    end
    cdeq = cd - gd * vd;
    A[n1p]-= cdeq
    A[n2] += cdeq
end
 

function nameOfMatrixEntries(device::Diode)
    return device.terminals
end
