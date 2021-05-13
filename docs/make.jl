push!(LOAD_PATH,"../src/")
using Documenter, DisastersDouenne

makedocs(modules = [DisastersDouenne], sitename = "DisastersDouenne.jl")

deploydocs(repo = "github.com/remihannotel/DisastersDouenne.jl.git", devbranch = "main")
