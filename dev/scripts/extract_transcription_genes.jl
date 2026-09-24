#!/usr/bin/env julia
#
# Generate `src/organisms/coreA/data/transcription_genes.tsv` — the per-gene
# sequence extract spec §11 task 10.1 vendors, one row for each of Core A′'s
# seventeen genes.
#
# WHY THIS FILE EXISTS. The published model rebuilds every transcription rate
# constant from a genome record, two spreadsheets and a counts table, none of
# which are vendored here (spec §5: raw upstream tables are not vendored, and
# compute nodes have no network). The rate constant needs only six numbers per
# gene — transcript length, four base counts and a protein copy number — plus
# the measured transcript mean the calibration check compares against. This
# reshapes those five sources into one tab-separated table and changes no value.
#
# WHAT LEGITIMATELY INVALIDATES THE OUTPUT. A different `Minimal_Cell` commit;
# a change to the seventeen loci Core A′ carries; a correction to the
# complement-strand handling below. Nothing else — in particular, re-running
# this script against the same checkout must reproduce the file byte for byte,
# which `extract_transcription_genes.slurm` checks.
#
# THE LOOKUP CHAIN, and why it takes five files rather than the three the
# archived design named (spec §12, amendment of 2026-09-10). Transcript length
# and the four base counts come from `syn3A.gb` alone. The protein copy number
# does not: `syn3A.gb` carries no AOE protein ids at all, so the count is
# reached the way `MinCell_CMEODE.py:57-110` and `:146-170` reach it —
#
# `!First2` carries the transcript's first two bases. The published rate law
# reads their NTP concentrations as C1 and C2 (`MinCell_CMEODE.py:370-372`,
# `CMono1`/`CMono2`), so the two bases are part of the per-gene data even though
# the archived design's header omitted them.
#
# `!Residues` is the protein length translation charges for (spec §11 task
# 11.1): the amino acids of the transcript translated under NCBI table 4, as
# `MinCell_CMEODE.py:121` builds `aasequence`, with the one stop excluded. It
# is counted from the translation rather than computed as `length/3 - 1`, and
# the two are asserted equal, so a CDS with an internal stop or a partial codon
# aborts instead of loading a wrong residue count (§12, 2026-09-23 D).
#
#   JCVISYN3A_xxxx  --suffix-->  MMSYN1_xxxx
#                   --Syn3A_annotation_compilation.xlsx col 6 -> col 14-->
#                   JCVSYN2_xxxxx
#                   --syn2.gb CDS /protein_id-->  AOE_xxxxx.x
#                   --proteomics.xlsx "Protein" -> col 22-->  copies
#
# Keying the annotation table on the MMSYN1 code is what the published model
# does; the JCVISYN3A and MMSYN1 tags share their numeric suffix. Every step
# fails loudly and by name, because a silently absent gene would show up only
# as a model with sixteen transcripts.
#
# Usage:
#   julia --project=dev/scripts dev/scripts/extract_transcription_genes.jl \
#       <Minimal_Cell checkout> <output .tsv>

using XLSX
using Printf
using SHA: sha256

# The revision this extract is derived from. The sibling extractors pin the
# same constant and refuse a mismatch, and for the same reason: regenerating
# against another revision would rewrite all seventeen genes' numbers, round-
# trip cleanly against the wrong source, and surface only as an unexplained
# `git diff`. The `.slurm` determinism check cannot catch it either, since it
# runs both passes against the same checkout.
const UPSTREAM_COMMIT = "db048aca5fe85438e0129819bbf0314b037dd931"

# Core A′'s seventeen genes: ten glycolytic, four phosphotransferase, three
# nucleotide-recycling (`dev/notes/reduced-syn3a-scoping.md`). The reaction
# names are carried here only to make an error message readable — the extract
# is keyed on the locus, which is the join key every consumer uses.
const GENES = [
    ("JCVISYN3A_0445", "PGI"), ("JCVISYN3A_0220", "PFK"),
    ("JCVISYN3A_0131", "FBA"), ("JCVISYN3A_0727", "TPI"),
    ("JCVISYN3A_0607", "GAPD"), ("JCVISYN3A_0606", "PGK"),
    ("JCVISYN3A_0729", "PGM"), ("JCVISYN3A_0213", "ENO"),
    ("JCVISYN3A_0221", "PYK"), ("JCVISYN3A_0475", "LDH_L"),
    ("JCVISYN3A_0233", "ptsI"), ("JCVISYN3A_0234", "Crr"),
    ("JCVISYN3A_0694", "ptsH"), ("JCVISYN3A_0779", "ptsG"),
    ("JCVISYN3A_0651", "ADK1"), ("JCVISYN3A_0344", "PPA"),
    ("JCVISYN3A_0203", "GK1"),
]

# `MinCell_CMEODE.py:85` — a protein with no quantification gets this rather
# than a zero, and the floor applies to quantified ones too via `max`.
const DEFAULT_PTN_COUNT = 10

struct CDS
    first::Int
    last::Int
    complement::Bool
    protein_id::Union{String, Nothing}
end

"""
    read_genbank(path) -> (Dict{String,CDS}, String)

The CDS features keyed by locus tag, and the ORIGIN sequence uppercased.

Only the two location forms this file actually uses are accepted — `a..b` and
`complement(a..b)`. A `join(...)` location would be silently mis-extracted, so
it raises instead; neither genome record contains one.
"""
function read_genbank(path::AbstractString)
    features = Dict{String, CDS}()
    seq = IOBuffer()
    in_origin = false
    in_cds = false
    seen_qualifier = false
    loc = ""
    tag = nothing
    pid = nothing

    flush_cds!() = begin
        if in_cds && tag !== nothing
            m = match(r"^(complement\()?(\d+)\.\.(\d+)\)?$", loc)
            m === nothing && error(
                "$(basename(path)): CDS $tag has location '$loc', which is " *
                "neither 'a..b' nor 'complement(a..b)'. This script extracts " *
                "only those two forms; a compound location would be " *
                "mis-extracted rather than reported")
            features[tag] = CDS(parse(Int, m[2]), parse(Int, m[3]),
                                m[1] !== nothing, pid)
        end
        in_cds = false; seen_qualifier = false; tag = nothing; pid = nothing; loc = ""
    end

    for line in eachline(path)
        if in_origin
            startswith(line, "//") && break
            for c in line
                (c in ('a', 'c', 'g', 't', 'A', 'C', 'G', 'T', 'n', 'N')) &&
                    write(seq, uppercase(c))
            end
            continue
        end
        if startswith(line, "ORIGIN")
            flush_cds!()
            in_origin = true
            continue
        end
        # A feature key sits at column 6; its qualifiers at column 22.
        if length(line) > 5 && line[1:5] == "     " && line[6] != ' '
            flush_cds!()
            parts = split(strip(line), r"\s+"; limit = 2)
            if parts[1] == "CDS"
                in_cds = true
                seen_qualifier = false
                loc = length(parts) > 1 ? parts[2] : ""
            end
        elseif in_cds && startswith(line, "                     ")
            q = strip(line)
            if startswith(q, "/")
                seen_qualifier = true
                if startswith(q, "/locus_tag=")
                    tag = strip(q[12:end], ['"'])
                elseif startswith(q, "/protein_id=")
                    pid = strip(q[13:end], ['"'])
                end
            elseif !seen_qualifier
                # A continuation of the location line, which happens only for
                # compound locations. Keep it so flush_cds! reports it. Once a
                # qualifier has started, an unprefixed line is that qualifier's
                # wrapped value (/translation is hundreds of characters) and
                # must not be mistaken for part of the location.
                loc *= q
            end
        end
    end
    flush_cds!()
    return features, String(take!(seq))
end

const COMPLEMENT = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'N' => 'N')

# NCBI translation table 4 (mycoplasma and spiroplasma), codons ordered with
# U, C, A, G varying first-base-slowest. It is the standard code with UGA read
# as tryptophan, which is why a syn3A gene cannot be translated with table 1.
const TABLE4_AA = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
const RNA_ORDER = Dict('U' => 0, 'C' => 1, 'A' => 2, 'G' => 3)

"""
    residues(rna) -> Int

The number of amino acids translation charges for: the transcript translated
under table 4, less its one stop codon. Raises on a partial codon, a base
outside ACGU, or any stop other than a single terminal one, since each would
make `length/3 - 1` and the translated count disagree.
"""
function residues(rna::AbstractString)
    n = length(rna)
    n % 3 == 0 || error("transcript of $n nucleotides is not a whole number of codons")
    aa = Char[]
    for i in 1:3:n
        c = rna[i:i+2]
        all(haskey(RNA_ORDER, b) for b in c) ||
            error("codon $c at position $i holds a base outside ACGU")
        push!(aa, TABLE4_AA[16 * RNA_ORDER[c[1]] + 4 * RNA_ORDER[c[2]] + RNA_ORDER[c[3]] + 1])
    end
    stops = findall(==('*'), aa)
    stops == [length(aa)] || error(
        "translation under table 4 has stops at codons $stops of $(length(aa)); " *
        "expected exactly one, at the end")
    r = length(aa) - 1
    r == n ÷ 3 - 1 || error("translated residues $r disagree with length/3 - 1")
    return r
end

"""
    transcript(cds, genome) -> String

The mature transcript, as the published model builds it: extract the feature's
span, reverse-complement it when the feature is on the complement strand, then
transcribe T to U.
"""
function transcript(c::CDS, genome::AbstractString)
    dna = genome[c.first:c.last]
    if c.complement
        dna = String([COMPLEMENT[ch] for ch in Iterators.reverse(dna)])
    end
    return replace(dna, 'T' => 'U')
end

"""
    annotation_map(path) -> Dict{String,String}

MMSYN1 locus tag to JCVSYN2 locus tag, from the condensed annotation sheet.
Columns 6 and 14 are `MinCell_CMEODE.py:587`'s `iloc[:,5]` and `iloc[0,13]`.
"""
function annotation_map(path::AbstractString)
    d = XLSX.getdata(XLSX.readxlsx(path)["Syn3A_annotation_compilation_condensed"])
    out = Dict{String, String}()
    for i in 2:size(d, 1)
        mm = d[i, 6]
        j2 = d[i, 14]
        (mm isa AbstractString && j2 isa AbstractString) || continue
        k = strip(mm)
        # `iloc[0]` in the published loader: the FIRST row for a key wins. A
        # Dict assignment would keep the last, so a duplicated key would hand
        # us a different value from the model we reproduce. Refuse instead.
        haskey(out, k) && out[k] != strip(j2) && error(
            "$(basename(path)): two rows for $k give $(out[k]) and " *
            "$(strip(j2)). The published loader takes the first; resolve the " *
            "duplicate rather than letting row order decide.")
        haskey(out, k) || (out[k] = strip(j2))
    end
    return out
end

"""
    proteomics_map(path) -> Dict{String,Float64}

AOE protein id to copy number. `pandas` reads this with `skiprows=[0]`, so the
header is spreadsheet row 2 and the data start at row 3; `iloc[0,21]` is
spreadsheet column 22.
"""
function proteomics_map(path::AbstractString)
    d = XLSX.getdata(XLSX.readxlsx(path)["Proteomics"])
    out = Dict{String, Float64}()
    for i in 3:size(d, 1)
        p = d[i, 1]
        v = d[i, 22]
        (p isa AbstractString && v isa Real) || continue
        k = strip(p)
        # Same first-wins rule as above; replicate measurements for one protein
        # are the obvious way this could bite.
        haskey(out, k) && out[k] != float(v) && error(
            "$(basename(path)): two rows for $k give $(out[k]) and $(float(v)). " *
            "The published loader takes the first; resolve the duplicate.")
        haskey(out, k) || (out[k] = float(v))
    end
    return out
end

"""
    mrna_map(path) -> Dict{String,Float64}

Syn3A locus tag to measured mean transcript count, from the plain CSV.
"""
function mrna_map(path::AbstractString)
    out = Dict{String, Float64}()
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    tag_col = findfirst(==("LocusTag"), header)
    cnt_col = findfirst(==("Count"), header)
    (tag_col === nothing || cnt_col === nothing) && error(
        "$(basename(path)): expected 'LocusTag' and 'Count' columns, got $header")
    for line in lines[2:end]
        isempty(strip(line)) && continue
        f = split(strip(line), ',')
        length(f) < max(tag_col, cnt_col) && continue
        k = strip(f[tag_col])
        v = parse(Float64, f[cnt_col])
        # Same first-wins rule as above.
        haskey(out, k) && out[k] != v && error(
            "$(basename(path)): two rows for $k give $(out[k]) and $v. The " *
            "published loader takes the first; resolve the duplicate.")
        haskey(out, k) || (out[k] = v)
    end
    return out
end

function main()
    length(ARGS) == 2 || error(
        "usage: extract_transcription_genes.jl <Minimal_Cell checkout> <output .tsv>")
    checkout, outpath = ARGS
    md = joinpath(checkout, "CME_ODE", "model_data")
    isdir(md) || error("no CME_ODE/model_data under '$checkout'")

    # A checkout that is not a git repository is warned about rather than
    # refused, since a tarball is a legitimate way to have the files; a
    # checkout at the wrong revision is refused outright.
    # The full hash, not `--short`: an abbreviation's length is a git config
    # setting, so `core.abbrev = 12` would otherwise produce a different file
    # from the same checkout and break the byte-for-byte round trip.
    head = try
        strip(read(pipeline(`git -C $checkout rev-parse HEAD`,
                            stderr = devnull), String))
    catch
        nothing
    end
    if head === nothing
        @warn "Could not read a git revision from $checkout, so the upstream " *
              "commit is unverified. Expected $UPSTREAM_COMMIT."
    elseif head != UPSTREAM_COMMIT
        error("$checkout is at $head, not the $UPSTREAM_COMMIT this extract " *
              "is derived from. Check the revision out, or update " *
              "UPSTREAM_COMMIT here and in src/organisms/coreA/data/README.md.")
    end
    commit = UPSTREAM_COMMIT

    syn3a, genome3a = read_genbank(joinpath(md, "syn3A.gb"))
    syn2, _ = read_genbank(joinpath(md, "syn2.gb"))
    ann = annotation_map(joinpath(md, "FBA", "Syn3A_annotation_compilation.xlsx"))
    prot = proteomics_map(joinpath(md, "proteomics.xlsx"))
    mrna = mrna_map(joinpath(md, "mRNA_counts.csv"))

    rows = String[]
    totals = Dict('A' => 0, 'C' => 0, 'G' => 0, 'U' => 0)

    for (locus, reaction) in GENES
        haskey(syn3a, locus) || error(
            "locus $locus ($reaction) is absent from syn3A.gb. A missing gene " *
            "would otherwise surface only as a model with $(length(GENES) - 1) " *
            "transcripts")
        rna = transcript(syn3a[locus], genome3a)
        n = length(rna)
        counts = Dict(b => count(==(b), rna) for b in ('A', 'C', 'G', 'U'))
        s = sum(values(counts))
        s == n || error(
            "locus $locus ($reaction): base counts sum to $s but the transcript " *
            "is $n nucleotides long, so it holds a base outside ACGU")
        for b in ('A', 'C', 'G', 'U')
            totals[b] += counts[b]
        end

        mm = "MMSYN1_" * split(locus, '_')[2]
        haskey(ann, mm) || error(
            "locus $locus ($reaction) maps to $mm, which is absent from " *
            "Syn3A_annotation_compilation.xlsx")
        j2 = ann[mm]
        haskey(syn2, j2) || error(
            "locus $locus ($reaction) maps to $mm -> $j2, which is absent " *
            "from syn2.gb")
        aoe = syn2[j2].protein_id
        aoe === nothing && error(
            "locus $locus ($reaction) maps to $j2 in syn2.gb, which carries " *
            "no /protein_id, so it has no AOE id to look up")
        haskey(prot, aoe) || error(
            "locus $locus ($reaction) maps to AOE id $aoe, which is absent " *
            "from proteomics.xlsx")
        ptn = max(DEFAULT_PTN_COUNT, round(Int, prot[aoe]))

        haskey(mrna, locus) || error(
            "locus $locus ($reaction) is absent from mRNA_counts.csv")

        push!(rows, join([locus, string(n), string(counts['A']),
                          string(counts['C']), string(counts['G']),
                          string(counts['U']), rna[1:2], string(ptn),
                          @sprintf("%.4f", mrna[locus]), string(residues(rna)),
                          "syn3A.gb|Syn3A_annotation_compilation.xlsx|syn2.gb|proteomics.xlsx|mRNA_counts.csv"],
                         '\t'))
    end

    mkpath(dirname(outpath))
    open(outpath, "w") do io
        println(io, "!!SBtab TableType='Quantity' TableName='transcription genes' Document='coreA'")
        println(io, "% Generated by dev/scripts/extract_transcription_genes.jl from")
        println(io, "% Luthey-Schulten-Lab/Minimal_Cell at commit $commit. Do not edit by hand.")
        println(io, "% Base-count totals across the seventeen genes: " *
                    "A $(totals['A']), C $(totals['C']), G $(totals['G']), U $(totals['U']).")
        println(io, join(["!ID", "!Length", "!A", "!C", "!G", "!U", "!First2",
                          "!PtnCount", "!MeanMRNA", "!Residues", "!UpstreamRow"], '\t'))
        for r in rows
            println(io, r)
        end
    end

    println("wrote $outpath: $(length(rows)) rows")
    println("base totals: A $(totals['A']), C $(totals['C']), " *
            "G $(totals['G']), U $(totals['U'])")
    # Printed so the README can record it: a hand-edit to a value no test
    # asserts is invisible otherwise.
    println("sha256: ", bytes2hex(sha256(read(outpath))))
end

main()
