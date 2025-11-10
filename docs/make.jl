using AxialAM
using Documenter

DocMeta.setdocmeta!(AxialAM, :DocTestSetup, :(using AxialAM); recursive=true)

makedocs(;
    modules=[AxialAM],
    authors="Fernando Peñaranda",
    sitename="AxialAM.jl",
    format=Documenter.HTML(;
        canonical="https://fernandopenaranda.github.io/AxialAM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fernandopenaranda/AxialAM.jl",
    devbranch="main",
)
