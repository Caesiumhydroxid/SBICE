# SBICE - Speedy* Bayesian Inference for Circuit Estimation  🎲 🔌⚡

SBICE is a probabilistic electrical circuit simulator that can be used to get (hopefully better) estimates of the 
actual parameters of each component in the circuit.

*results may vary ([SPICE](https://de.wikipedia.org/wiki/SPICE_(Software)))

# Supported Components

- Resistors
- Capacitors
- Inductors
- Voltage Sources 
- Current Sources
- Diodes (only DC Model) 
- Voltage controlled voltage sources (as a simple model for operational amplifiers)

# Setting up inference 🧑‍💻

## Required information🗄️
- SPICE-like input file of the circuit ([example](measurements/sallenkeyExt/sallenkeyx.cir))
- Tolerances File ([example](measurements/sallenkeyExt/sallenkeyx.tol))
- Measurement Data for some nodes in the circuit (either DC-Voltage-Sweep or AC-frequency-sweep) ([example](measurements/sallenkeyExt) SALEX01-SALEX.CSV)

## Run the analysis

In the folder [analyses](analyses) there are several files which were used while developing this project which should help getting started.
- Sampling versions
    - [Linear DC Analysis](analyses/divider.jl)
    - [AC Analysis](analyses/salkxLess.jl)
    - [Non-Linear DC Analysis](analyses/diodean.jl)
- VI versions
  - [Linear DC Analysis](analyses/dividerADVI.jl)
  - [AC Analysis](analyses/salkxLessADVI.jl)
  - [Non-Linear DC Analysis](analyses/diodeADVI.jl)

## Adjusting the Priors 
The probabilistic models are located in [this file](probabilistic.jl) which contain all the chosen distributions for *device parameters* and measurement *noise* and *offsets*. This set of choices might not be optimal for any circuit but it should be clear how to modify these.

# Results

The results are somewhat promising depending on the quality of the measurements. 

# TODO's
- Implement more devices
- Make it easier to use (additional utility functions to create the input to the model)
- Try out different samplers (especially for the fault model with discrete choice variables)
- Make the non-linear solver numerically stable 

# Remarks
This project was created for the course 194.150 Probabilistic Programming and AI at TU Wien.
