using Documenter
using InferCell

DocMeta.setdocmeta!(InferCell, :DocTestSetup, :(using InferCell); recursive = true)

makedocs(;
    modules = [InferCell],
    authors = "Tom Kimpson",
    sitename = "InferCell.jl",
    format = Documenter.HTML(;
        canonical = "https://macsys-coe.github.io/InferCell",
        edit_link = "main",
        repolink = "https://github.com/MACSYS-CoE/InferCell",
        assets = String[],
        # The prose carries rate laws and conservation residuals; KaTeX ships
        # with Documenter, which is why docs/javascripts/mathjax.js went away.
        mathengine = Documenter.KaTeX(),
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "The reduced syn3A model" => "user-guide/reduced-syn3a.md",
        "Toy models" => [
            "ODE inference (NUTS)" => "user-guide/ode-inference.md",
            "SSA inference (ABC-SMC)" => "user-guide/ssa-inference.md",
            "Boundary protocol" => "user-guide/boundary-protocol.md",
        ],
        "API reference" => [
            "Overview" => "api/index.md",
            "Parameters" => "api/parameters.md",
            "Sub-model interface" => "api/interface.md",
            "Interface contract" => "api/corea-interface.md",
            # One page overflowed Documenter's 200 KiB HTML limit, so the
            # modules' docstrings are split by block (KIM-22).
            "Model modules" => [
                "Overview" => "api/corea-modules.md",
                "Metabolism" => "api/corea-metabolism.md",
                "Gene expression" => "api/corea-gene-expression.md",
                "Assembly, priors and checks" => "api/corea-assembly.md",
            ],
            "Observables and data" => "api/corea-observables.md",
            "Orchestrator" => "api/orchestrator.md",
            "Handshake driver" => "api/handshake.md",
            "Conditional sampler" => "api/conditional-sampler.md",
            "Inference (toy models)" => "api/inference.md",
            "Boundary (toy models)" => "api/boundary.md",
            "Toy models" => "api/models.md",
        ],
    ],
    # An exported symbol that no `@docs` block includes fails the build. This
    # is the whole point of the migration: the hand-written API pages drifted
    # from the code (a demoted check was still documented as a throw), and a
    # signature can no longer be copied by hand and left behind.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/MACSYS-CoE/InferCell",
    devbranch = "main",
)
