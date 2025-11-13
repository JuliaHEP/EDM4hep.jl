using Documenter
using EDM4hep
using EDM4hep.RootIO
using Literate

name = ARGS[1]

examplesdir =  joinpath(@__DIR__, "src/examples")
file = "$(name).lit"

println("====> Processing $name")
if haskey(ENV, "GITHUB_ACTIONS")
    println("Running inside GitHub Actions CI, skipping execution of notebooks")
else
    Literate.notebook(joinpath(examplesdir, file), examplesdir, name = name, execute = false, documenter = true, credit = true)
    Literate.script(joinpath(examplesdir, file), examplesdir, name = name, keep_comments = false, documenter = true, credit = false)
    Literate.markdown(joinpath(examplesdir, file), examplesdir, name = name, execute = true, documenter = true, continue_on_error=true, credit = true)
end

