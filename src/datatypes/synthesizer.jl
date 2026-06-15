@data BenchmarkType begin
    None
    UntilFound
    UntilTimeout
end

@data Iterator begin
    BreadthFirst
    DepthFirst
    MaxBond
    DeltaEnergy
    BottomUp
end

struct SynthesizerSettings
    max_programs::Int
    max_time::Int
    max_depth::Int
    iterator::Iterator
    benchmark_type::BenchmarkType
    goal::Dict{Any, Bool}
    options::Dict{Symbol, Any}
end
function SynthesizerSettings(;
        max_programs::Int = typemax(Int),
        max_time::Int = typemax(Int),
        max_depth::Int = typemax(Int),
        iterator::Iterator = BreadthFirst,
        benchmark_type::BenchmarkType = None,
        goal::Any = Dict{Any, Bool}(),
        options::Dict{Symbol, Any} = Dict{Symbol, Any}()
)
    if !(goal isa Dict{Any, Bool})
        new_goal = Dict{Any, Bool}()
        if goal isa Vector
            for item in goal
                new_goal[item] = false
            end
        else
            throw("Goal must be a Dict{Any, Bool} or a Vector.")
        end
        goal = new_goal
    end
    return SynthesizerSettings(
        max_programs, max_time, max_depth, iterator, benchmark_type, goal, options
    )
end

function get_iterator(
        settings::SynthesizerSettings, grammar::ContextSensitiveGrammar, starting_symbol::Symbol
)
    @match settings.iterator begin
        # BreadthFirst => return BFSIterator(grammar, starting_symbol, max_depth=settings.max_depth)
        BreadthFirst => return MyBFSIterator(
            grammar, starting_symbol, max_depth = settings.max_depth)
        DepthFirst => return DFSIterator(
            grammar, starting_symbol, max_depth = settings.max_depth)
        BottomUp => return BUIterator(
            grammar, starting_symbol, max_depth = settings.max_depth)
        MaxBond => return MaxBondIterator(
            grammar, starting_symbol, max_depth = settings.max_depth)
        DeltaEnergy => return DeltaEnergyIterator(
            grammar, starting_symbol, max_depth = settings.max_depth)
    end
end

function check_stop_condition(
        settings::SynthesizerSettings,
        start_time::Float64,
        all_candidates,
        last_candidate;
        check_all_candidates = false
)
    if length(all_candidates) >= settings.max_programs
        return true
    end
    if time() - start_time > settings.max_time
        return true
    end
    if settings.benchmark_type == UntilFound
        if check_all_candidates
            for item in all_candidates
                if haskey(settings.goal, item)
                    settings.goal[item] = true
                end
            end
        elseif haskey(settings.goal, last_candidate)
            settings.goal[last_candidate] = true
        end

        if all(values(settings.goal))
            return true
        end
    end
    return false
end

function find_programs!(
        iterator::ProgramIterator,
        settings::SynthesizerSettings,
        interpreter,
        candidates;
        start_time = time(),
        job = nothing
)
    check_stop_condition(
        settings, start_time, candidates, nothing; check_all_candidates = true
    )
    for (i, program) in enumerate(iterator)
        result = interpreter(program)
        push!(candidates, result)

        if !isnothing(job)
            update!(job)
            yield()
        end

        if check_stop_condition(settings, start_time, candidates, program)
            break
        end
    end

    return iterator
end
