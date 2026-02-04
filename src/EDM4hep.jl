"""
Main module for `EDM4hep.jl` -- Key4hep Event Data Model for Julia.

All data model types are exported from this module for public use

# Exports

"""
module EDM4hep
    using Preferences

    const edmodel = get(ENV, "EDModel", Preferences.load_preference(@__MODULE__, "edmodel", "hep"))

    @show edmodel  

    function set_edmodel(new_edmodel::String)
        if !(new_edmodel in ("hep", "eic"))
            throw(ArgumentError("Invalid edmodel: \"$(new_edmodel)\""))
        end
        # Set it in our runtime values, as well as saving it to disk
        @set_preferences!("edmodel" => new_edmodel)
        @info("New edmodel set; restart your Julia session for this change to take effect!")
    end
    function get_edmodel()::String
        return edmodel
    end
    
    include("Components.jl")
    include("Datatypes.jl")
    include("EDCollection.jl")
    include("RootIO.jl")
    include("SystemOfUnits.jl")
    include("Histograms.jl")
    include("Analysis.jl")
    
end # module EDM4hep