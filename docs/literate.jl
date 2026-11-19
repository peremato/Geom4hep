#  Convert a Literate example (docs/src/examples/<name>.lit) into a Markdown page for Documenter
#  (executed, with its output and plots), a Jupyter notebook and a plain Julia script.
#  Usage: julia --project=docs docs/literate.jl <name>
using Literate
name = ARGS[1]

examplesdir = joinpath(@__DIR__, "src/examples")
println("====> Processing $name")
file = joinpath(examplesdir, "$(name).lit")

Literate.markdown(file, examplesdir, name = name, execute = true, documenter = true, credit = true)
Literate.notebook(file, examplesdir, name = name, execute = false, documenter = true, credit = true)
Literate.script(file, examplesdir, name = name, keep_comments = false, documenter = true, credit = false)
