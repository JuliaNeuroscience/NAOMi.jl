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
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "API" => [
            "Parameters" => "parameters.md",
            "Time traces" => "timetraces.md",
            "Optics" => "optics.md",
            "Volume" => "volume.md",
            "Scanning" => "scanning.md",
            "I/O" => "io.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/JuliaNeuroscience/NAOMi.jl",
    devbranch="main",
)
