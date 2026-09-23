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

# ---------------------------------------------------------------------------
# The published rate law's constants: `translation_rate_restart.py`, the law
# the 60 s rebuild runs for every minute of the cycle after the first.
# ---------------------------------------------------------------------------

"""
Ribosome turnover, `riboKcat`, in residues per second.

12, as both translation-rate files that run set it
(`translation_rate_start.py:18`, `translation_rate_restart.py:19`). The 10 at
`MinCell_CMEODE.py:332` is a module variable no rate law reads: `TranslatRate`
is imported from the rate files, which define their own.
"""
const RIBO_KCAT = 12.0

"`riboK0`, mM (`translation_rate_restart.py:20`)."
const RIBO_K0 = 4 * 25e-6

"""
`riboKd`, mM: 1e-3, the restart file's value (`translation_rate_restart.py:21`).

The start file has 1e-4, and it governs only the first minute. The published
model runs `MinCell_CMEODE.py` for one minute and then `MinCell_restart.py`
for every minute after, and the restart file's law is the one the 60 s
rebuild evaluates. This module's rebuild is that channel, so it takes the
restart law whole: this `riboKd` and the restart form of `kcat_mod` in
[`translation_kcat`](@ref). Chosen 2026-09-24 and recorded in spec §12.
"""
const RIBO_KD = 1e-3

"Ribosomes per cell, converted to mM at the registry's initial volume."
const RIBOSOME_COPIES = 503

"""
The number of per-amino-acid charged-tRNA pools the lumped pool stands in for.

The published law sums `n_aa · riboKd / [aa-tRNA]` over twenty pools and adds
`riboKd² / [fMet-tRNA]²`. Core A′ has one lumped charged pool, and each of
those twenty-one concentrations reads its per-amino-acid share,
`[M_trna_chg_c] / 20`. That keeps the published law's structure: at the
nominal 0.2 mM charged pool the share is 0.01 mM against upstream's 150
copies, 0.0074 mM. Substituting the whole pool for each per-amino-acid pool
would move each term's concentration 27-fold, and with it the pool's
elasticity, from about 0.09 to about 0.005, for no reason but the lumping.
**Ours**, chosen 2026-09-24 (spec §12), and `reduction_notes` says so.
"""
const TL_AA_TYPES = 20

"Ribosomes per transcript never exceed this (`translation_rate_restart.py:103`)."
const POLYSOME_CAP = 15

"""
    ribosomes_per_transcript(length) -> Int

The published polysome size, `min(15, max(1, round(length/125 − 1)))`.
`round` ties to even in both Python 3 and Julia, and no transcript length
lands on a tie, since that would need `length/125` to end in .5.
"""
ribosomes_per_transcript(n::Integer) = min(POLYSOME_CAP, max(1, round(Int, n / 125 - 1)))

"""
    translation_kcat(length) -> Float64

The restart law's `kcat_mod`: `(0.25·n + 0.2)·riboKcat` for a polysome of
`n > 1` ribosomes, and `0.45·riboKcat` for one. The start file's `+0.25` would
apply for the first minute only (see [`RIBO_KD`](@ref)).
"""
function translation_kcat(n::Integer)
    r = ribosomes_per_transcript(n)
    return r > 1 ? (0.25 * r + 0.2) * RIBO_KCAT : 0.45 * RIBO_KCAT
end

"""
    translation_rate_constant(length, residues, chg_mM; ribo_conc) -> Float64

The published restart law (`translation_rate_restart.py:46-115`), equation 3,
with the lumped pool's per-amino-acid share `c = chg_mM / 20` standing in for
each of the twenty-one charged-tRNA concentrations:

```
k = kcat_mod / ( (1 + K₀/[ribosome])·K_d²/c² + residues·K_d/c + residues )
```

`residues` is the published `n_tot − 1`, since `n_tot` counts the stop codon,
and it is also `Σ n_aa`, since the stop is no amino acid. So the per-gene
composition drops out of the lumped law, which needs only the length and the
residue count.

The share is floored at one particle, at the registry's volume, as upstream
floors each pool with `max(1, count)`. Without it an exhausted charged pool
would divide by zero. `ribo_conc` is frozen at the registry's volume, as
upstream computes `ribosomeConc` once at the initial radius.

Written once, as a free function, so tests assert the module's own arithmetic.
"""
function translation_rate_constant(n::Integer, residues::Integer, chg_mM::Real;
                                   ribo_conc::Real = RIBOSOME_COPIES /
                                                     corea_particles_per_mM())
    c = max(chg_mM / TL_AA_TYPES, 1 / corea_particles_per_mM())
    denom = (1 + RIBO_K0 / ribo_conc) * RIBO_KD^2 / c^2 +
            residues * RIBO_KD / c + residues
    return translation_kcat(n) / denom
end

"""
    translation_elasticity(length, residues, chg_mM; ribo_conc) -> Float64

`∂ ln k / ∂ ln [M_trna_chg_c]`, analytically: the charged-pool terms' share of
the denominator, with the fMet term counted twice because it is squared. This
is the charged-tRNA reverse channel's gain (spec §9, R1, task 11.10).
"""
function translation_elasticity(n::Integer, residues::Integer, chg_mM::Real;
                                ribo_conc::Real = RIBOSOME_COPIES /
                                                  corea_particles_per_mM())
    c = chg_mM / TL_AA_TYPES
    a = (1 + RIBO_K0 / ribo_conc) * RIBO_KD^2 / c^2
    b = residues * RIBO_KD / c
    return (2a + b) / (a + b + residues)
end

export read_translation_residues, ribosomes_per_transcript, translation_kcat,
       translation_rate_constant, translation_elasticity
export RIBO_KCAT, RIBO_K0, RIBO_KD, RIBOSOME_COPIES, TL_AA_TYPES, POLYSOME_CAP
