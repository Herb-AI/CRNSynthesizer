# CRNSynthesizer

Part of the [Research Project](https://github.com/TU-Delft-CSE/Research-Project) 2026 Q4 at [TU Delft](https://github.com/TU-Delft-CSE).

## Purpose
Automated discovery of chemical reaction networks based on (partial) measurements.

## Usage

### Installing the Package
To install the CRNSynthesizer package, open Julia in the project root directory and enter the package manager by pressing `]`, then run:

```julia
dev .
```

### Running Benchmarks
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


## Looking for more benchmarks
Suggestions for additional reaction networks to expand the benchmark suite are very welcome! Please open an issue or pull request if you have ideas or datasets to contribute.
