module SynRXNLoader

using PythonCall
using DataFrames

export load_synrxn_dataset, get_available_datasets

const _SYNRXN_CACHE_DIR = joinpath(homedir(), ".cache", "synrxn")

function _synrxn_dataloader(task::String; version = "v1.0.0", source = "github",
        gh_enable = true, cache_dir = _SYNRXN_CACHE_DIR, kwargs...)
    @pyconst(pyimport("synrxn.data").DataLoader)(
        task; version, source, gh_enable, cache_dir, kwargs...)
end

"""
    get_available_datasets(task; version="v1.0.0", source="github", gh_enable=true, kwargs...) -> Vector{String}

Return available dataset names for the given task family.
"""
function get_available_datasets(
        task::String; version = "v1.0.0", source = "github", gh_enable = true, kwargs...)
    pyconvert(
        Vector{String}, _synrxn_dataloader(
            task; version, source, gh_enable, kwargs...).available_names())
end

"""
    load_synrxn_dataset(task, name; version="v1.0.0", source="github", gh_enable=true, kwargs...) -> DataFrame

Load a SynRXN dataset as a Julia `DataFrame`.
"""
function load_synrxn_dataset(task::String, name::String; version = "v1.0.0",
        source = "github", gh_enable = true, kwargs...)
    py_df = _synrxn_dataloader(task; version, source, gh_enable, kwargs...).load(name)
    return DataFrame(PyTable(py_df))
end

end
