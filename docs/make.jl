using NAOMi
using Documenter

DocMeta.setdocmeta!(NAOMi, :DocTestSetup, :(using NAOMi); recursive=true)

makedocs(;
    modules=[NAOMi],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    sitename="NAOMi.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaNeuroscience.github.io/NAOMi.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaNeuroscience/NAOMi.jl",
    devbranch="main",
)
