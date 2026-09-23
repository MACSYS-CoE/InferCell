"""
Core A′ translation: seventeen reactions, one per transcript, plus the
translocation of ptsG into the membrane. Spec §11 phase 11.

The second Core A′ module in the stochastic block. It reads the transcripts
[`CoreATranscription`](@ref) owns and publishes the protein counts the rest of
the model reads: the thirteen enzyme counts fill the metabolic modules' enzyme
slots through catalytic edges those modules declare (phase 11a), and the four
phosphotransferase carriers are credited to [`PtsTransport`](@ref)'s states
through producer counters.
"""

"""
    read_translation_residues([path]) -> Dict{Symbol,Int}

Each locus's residue count, from the `!Residues` column of the gene extract
phase 10 vendors (`src/organisms/coreA/data/transcription_genes.tsv`).

A residue is an amino acid translation charges for: the transcript translated
under NCBI table 4, less its one stop codon. The extract's generator counts the
translation and asserts it equals `length/3 − 1`, so this reader asserts the
same rather than trusting either (spec §12, 2026-09-23 D). The column is read
by name, not position, so phase 10's positional columns are untouched.
"""
function read_translation_residues(path::AbstractString = TRANSCRIPTION_EXTRACT)
    header = nothing
    out = Dict{Symbol, Int}()
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "!!") || startswith(s, "%")) && continue
        f = split(s, '\t')
        if startswith(s, "!")
            header = [strip(c, '!') for c in f]
            continue
        end
        col = header === nothing ? nothing : findfirst(==("Residues"), header)
        col === nothing && throw(ArgumentError(
            "$path has no !Residues column. Regenerate it with " *
            "dev/scripts/extract_transcription_genes.jl; see " *
            "src/organisms/coreA/data/README.md"))
        out[Symbol(f[1])] = parse(Int, f[col])
    end
    for g in read_transcription_genes(path)
        haskey(out, g.locus) || throw(ArgumentError(
            "$path holds no residue count for $(g.locus) ($(g.reaction))"))
        out[g.locus] == g.length ÷ 3 - 1 && g.length % 3 == 0 || throw(ArgumentError(
            "$(g.locus) ($(g.reaction)): $(out[g.locus]) residues, but a " *
            "$(g.length)-nucleotide transcript with one terminal stop has " *
            "$(g.length ÷ 3 - 1)"))
    end
    return out
end

export read_translation_residues
