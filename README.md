# CRNSynthesizer (Energy-Guided Search Extension)

This code is developed as part of the 2026 Q4 [Research Project](https://github.com/TU-Delft-CSE/Research-Project) of [TU Delft](https://github.com/TU-Delft-CSE).

## About this extension

This project extends [CRNSynthesizer](https://github.com/RichardWijers/CRNSynthesizer) by Richard Wijers.

This extension adds two new search heuristics that use bond-breaking energy to guide the search:

- **Maximum bond-order heuristic** - prioritises reactions with lower maximum bond order in reactants 
- **Delta energy heuristic** - prioritises reactions with lower net energy change

Both heuristics are benchmarked against baseline breadth-first search across seven CRNs (water formation, methane combustion, photosynthesis, ethylene glycol formation, methyl acetate hydrolysis, esterification, and glucose fermentation).

The full set of experiments comparing these heuristics against BFS is implemented in [`benchmark/heuristics_comparison.jl`](benchmark/heuristics_comparison.jl).


### Usage

The new heuristics are used as drop-in replacements for the default iterator in `SynthesizerSettings`. Specify `MaxBond` or `DeltaEnergy` as the `iterator`. For example:

```julia
synthesis_settings = SynthesizerSettings(;goal = [network_goal], iterator = MaxBond)
```

See [`benchmark/heuristics_comparison.jl`](benchmark/heuristics_comparison.jl) to see it be implemented

---

This work is an extension of [CRNSynthesizer](https://github.com/RichardWijers/CRNSynthesizer) by Richard Wijers. The original README is reproduced below:


## CRNSynthesizer

### Purpose
Automated discovery of chemical reaction networks based on (partial) measurements.

### Usage

#### Installing the Package
To install the CRNSynthesizer package, open Julia in the project root directory and enter the package manager by pressing `]`, then run:

```julia
dev .
```

#### Running Benchmarks
To run the provided benchmarks, use Julia from the command line:

```sh
julia benchmark/accuracy.jl
julia benchmark/feasibility.jl
```

<!-- ### Running a Search for Your Own Measurements
To use CRNSynthesizer with your own (partial) measurements:

1. Prepare your measurement data in a format similar to the examples in `benchmark/data/`.
2. Write a Julia script that loads your data and runs the synthesizer. For example:

```julia
using CRNSynthesizer
# Load your data (replace with your file)
data = include("benchmark/data/your_data.jl")
# Run the synthesizer
result = CRNSynthesizer.synthesize(data)
println(result)
```

3. Run your script with Julia:

```sh
julia your_script.jl
``` -->


### Looking for more benchmarks
Suggestions for additional reaction networks to expand the benchmark suite are very welcome! Please open an issue or pull request if you have ideas or datasets to contribute.
