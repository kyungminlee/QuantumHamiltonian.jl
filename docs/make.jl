using Documenter
using QuantumHamiltonian

makedocs(
    modules=[QuantumHamiltonian],
    doctest=true,
    sitename="QuantumHamiltonian.jl",
    format=Documenter.HTML(prettyurls=!("local" in ARGS)),
    authors="Kyungmin Lee",
    checkdocs=:all,
    pages = [
        "Home" => "index.md",
        "Types" => [
            "Hilbert space" => "hilbertspace.md",
            "Operator" => "operator.md",
            "Representation" => "representation.md",
            "Symmetry" => "symmetry.md"
        ],
        "Examples" => [
            "examples/spinhalf.md",
        ],
        "Index" => "links.md",
        "API" => "api.md",
    ]
  )

deploydocs(;
    repo="github.com/kyungminlee/QuantumHamiltonian.jl.git",
    devbranch="dev",
)
