module CRNSynthesizer

using HerbGrammar, HerbSpecification, HerbSearch, HerbInterpret, HerbConstraints, HerbCore
using MLStyle
using Catalyst, DiffEqParamEstim, Optimization, OptimizationNLopt, OrdinaryDiffEq
import RDKitMinimalLib, MoleculeFlow
using Term.Progress
using Term: with, update!

include("utils.jl")

include("datatypes/molecule.jl")
include("datatypes/reaction_network.jl")
include("datatypes/problem_definition.jl")
include("datatypes/grammar_data.jl")
include("datatypes/synthesizer.jl")

include("constraints/valid_molecule/uniform_ringbond.jl")
include("constraints/valid_molecule/generic_atom.jl")
include("constraints/valid_molecule/uniform_atom.jl")
include("constraints/valid_molecule/valid_SMILES.jl")
include("constraints/balanced_reaction/generic_balanced_reaction.jl")
include("constraints/balanced_reaction/uniform_balanced_reaction.jl")
include("constraints/balanced_reaction/balanced_reaction.jl")
include("constraints/contains_molecule.jl")
include("constraints/contains_reactions.jl")

include("grammars/bu_molecule_grammar.jl")
include("grammars/molecule_grammar.jl")
include("grammars/reaction_grammar.jl")
include("grammars/network_grammar.jl")

include("simulators/catalyst.jl")

include("interpreters/bu_molecule_interpreter.jl")
include("interpreters/molecule_interpreter.jl")
include("interpreters/reaction_interpreter.jl")
include("interpreters/network_interpreter.jl")

include("iterators/top_down.jl")
#include("iterators/bottom_up.jl")

include("synthesizers/similarity.jl")
include("synthesizers/molecule_synthesizer.jl")
include("synthesizers/reaction_synthesizer.jl")
include("synthesizers/network_synthesizer.jl")
include("synthesizers/fragment_composer.jl")

for n in names(@__MODULE__; all = true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end # module CRNSynthesizer
