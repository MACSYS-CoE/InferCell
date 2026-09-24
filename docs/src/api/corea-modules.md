```@meta
CurrentModule = InferCell
```

# Core A′ modules

The seven sub-models that make up the reduced syn3A organism. Each declares its
own coupling edges against the [Core A′ interface contract](corea-interface.md)
and vendors the parameter extract it imports from, so a module can be read,
tested and composed on its own.

Each module's vendored extracts live under `src/organisms/coreA/data/`, whose
`README.md` records the upstream file, the commit it was taken at and the
command that regenerates it byte for byte.

```@autodocs
Modules = [InferCell]
Pages = [
    "organisms/coreA/central_glycolysis.jl",
    "organisms/coreA/pts_transport.jl",
    "organisms/coreA/nucleotide_recycling.jl",
    "organisms/coreA/trna_charging.jl",
    "organisms/coreA/transcription.jl",
    "organisms/coreA/transcript_decay.jl",
    "organisms/coreA/translation.jl",
    "organisms/coreA/assembly.jl",
]
```
