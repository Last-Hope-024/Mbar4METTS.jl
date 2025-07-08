using Documenter

makedocs(
    sitename="Mbar4METTS Documentation",
    modules = [YourPackageName],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "https://github.com/Last-Hope-024/Mbar4METTS.jl.git",
)
