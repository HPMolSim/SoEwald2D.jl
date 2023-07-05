using SoEwald2D
using Documenter

DocMeta.setdocmeta!(SoEwald2D, :DocTestSetup, :(using SoEwald2D); recursive=true)

makedocs(;
    modules=[SoEwald2D],
    authors="Xuanzhao Gao",
    repo="https://github.com/ArrogantGao/SoEwald2D.jl/blob/{commit}{path}#{line}",
    sitename="SoEwald2D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ArrogantGao.github.io/SoEwald2D.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ArrogantGao/SoEwald2D.jl",
    devbranch="main",
)
