using NyquistCircuitModeling
using Documenter

DocMeta.setdocmeta!(NyquistCircuitModeling, :DocTestSetup, :(using NyquistCircuitModeling); recursive=true)

makedocs(;
    modules=[NyquistCircuitModeling],
    authors="Ir1n-a <159860234+Ir1n-a@users.noreply.github.com> and contributors",
    sitename="NyquistCircuitModeling.jl",
    format=Documenter.HTML(;
        canonical="https://Ir1n-a.github.io/NyquistCircuitModeling.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Ir1n-a/NyquistCircuitModeling.jl",
    devbranch="master",
)
