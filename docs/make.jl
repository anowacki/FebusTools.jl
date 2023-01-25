using FebusTools
using Documenter

DocMeta.setdocmeta!(FebusTools, :DocTestSetup, :(using FebusTools); recursive=true)

makedocs(;
    modules=[FebusTools],
    authors="Andy Nowacki <a.nowacki@leeds.ac.uk> and contributors",
    repo="https://github.com/anowacki/FebusTools.jl/blob/{commit}{path}#{line}",
    sitename="FebusTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://anowacki.github.io/FebusTools.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/anowacki/FebusTools.jl",
    devbranch="main",
)
