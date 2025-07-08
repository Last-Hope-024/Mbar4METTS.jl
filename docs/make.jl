using Documenter
using Mbar4METTS

makedocs(
    sitename = "Mbar4METTS",
    format = Documenter.HTML(),
    modules = [Mbar4METTS]
)


deploydocs(
    repo = "github.com/Last-Hope-024/Mbar4METTS.jl"
    ) 
