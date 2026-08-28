using Test
using InferCell
using Distributions

const FIXTURES = joinpath(@__DIR__, "fixtures")

@testset "Provenance-carrying parameter loader" begin

    central = read_source_table(joinpath(FIXTURES, "central_balanced.tsv");
                                file = "central_balanced")
    nucleotide = read_source_table(joinpath(FIXTURES, "nucleotide_balanced.tsv");
                                   file = "nucleotide_balanced")
    both = [central, nucleotide]

    @testset "Tables parse, skipping SBtab preamble and comments" begin
        @test central.file == "central_balanced"
        @test central.values["kcat_fwd_PGK3"] == 319.5
        @test central.values["kcat_fwd_ENO"] == 62.18
        @test nucleotide.values["kcat_fwd_PGK3"] == 140.8
        @test !haskey(central.values, "%")           # the comment line is skipped
        @test !haskey(central.values, "!ID")         # nor is the header
    end

    @testset "A missing file or column fails with a usable message" begin
        @test_throws ArgumentError read_source_table(joinpath(FIXTURES, "nope.tsv"))
        err = caught() do
            read_source_table(joinpath(FIXTURES, "central_balanced.tsv");
                              value_column = "NotAColumn")
        end
        @test err isa ArgumentError
        @test occursin("NotAColumn", err.msg)
        @test occursin("Mode", err.msg)      # lists the columns it does have
    end

    @testset "Informedness is derived from the uncertainty column" begin
        # A real balancing distribution.
        @test central.informedness["kcat_fwd_ENO"] == :balanced
        # The prior median at prior width: nothing informed it.
        @test central.informedness["conc_M_gtp_c"] == :prior_default
        # No quantified uncertainty at all, so any prior on it is ours.
        @test central.informedness["k_pts_GLCpts0_fwd"] == :asserted
    end

    @testset "Every loaded parameter names its file and identifier" begin
        p = load_parameter(both, "kcat_fwd_ENO";
                           name = :kcat_ENO, module_id = :Central,
                           prior = LogNormal(0.0, 1.0))
        @test p.value == 62.18
        @test source_file(p) == "central_balanced"
        @test provenance_of(p).identifier == "kcat_fwd_ENO"
        @test informedness(p) == :balanced
        @test !was_chosen_over_alternative(provenance_of(p))
    end

    @testset "An identifier in more than one file is reported" begin
        report = ambiguity_report(both)
        ids = [a.identifier for a in report]
        @test "kcat_fwd_PGK3" in ids
        @test "kcat_fwd_PYK3" in ids
        @test "conc_M_gtp_c" in ids
        # Present in exactly one file, so not ambiguous.
        @test "kcat_fwd_ENO" ∉ ids
        @test "kcat_fwd_ADK1" ∉ ids
    end

    @testset "Disagreements are distinguished from agreements" begin
        disagree = disagreements(both)
        ids = [a.identifier for a in disagree]

        # The kinetics trap and the concentration trap, the two errors of record.
        @test "kcat_fwd_PGK3" in ids
        @test "kcat_fwd_PYK3" in ids
        @test "conc_M_gtp_c" in ids
        @test length(ids) == 3

        # Same identifier, same value in both files: ambiguous but not a
        # disagreement.
        @test "kcat_fwd_SHARED" ∉ ids
        shared = only(filter(a -> a.identifier == "kcat_fwd_SHARED", ambiguity_report(both)))
        @test shared.agrees

        pgk3 = only(filter(a -> a.identifier == "kcat_fwd_PGK3", ambiguity_report(both)))
        @test !pgk3.agrees
        @test Set(pgk3.values) == Set(["central_balanced" => 319.5,
                                       "nucleotide_balanced" => 140.8])
    end

    @testset "The report is empty for a single-file load" begin
        # Emptiness is a positive result, not a failure to run.
        @test isempty(ambiguity_report([central]))
        @test isempty(disagreements([central]))
    end

    @testset "An ambiguous import without a governing file fails at load" begin
        err = caught() do
            load_parameter(both, "kcat_fwd_PGK3";
                           name = :kcat_PGK3, module_id = :Nucleotide,
                           prior = LogNormal(0.0, 1.0))
        end
        @test err isa ArgumentError
        # Names the identifier, both files and both values.
        @test occursin("kcat_fwd_PGK3", err.msg)
        @test occursin("central_balanced", err.msg)
        @test occursin("nucleotide_balanced", err.msg)
        @test occursin("319.5", err.msg)
        @test occursin("140.8", err.msg)
    end

    @testset "A declared governing file resolves the ambiguity" begin
        # PGK3 and PYK3 are nucleotide-module reactions, so the nucleotide file
        # governs — the correction of record.
        p = load_parameter(both, "kcat_fwd_PGK3";
                           name = :kcat_PGK3, module_id = :Nucleotide,
                           prior = LogNormal(0.0, 1.0),
                           governing = "nucleotide_balanced")
        @test p.value == 140.8
        @test source_file(p) == "nucleotide_balanced"
        @test was_chosen_over_alternative(provenance_of(p))
        @test provenance_of(p).alternatives == ["central_balanced" => 319.5]
    end

    @testset "A governing file that does not hold the identifier is rejected" begin
        err = caught() do
            load_parameter(both, "kcat_fwd_PGK3";
                           name = :kcat_PGK3, module_id = :Nucleotide,
                           prior = LogNormal(0.0, 1.0),
                           governing = "lipid_balanced")
        end
        @test err isa ArgumentError
        @test occursin("lipid_balanced", err.msg)
    end

    @testset "An identifier in no table is rejected" begin
        @test_throws ArgumentError load_parameter(both, "kcat_fwd_NOWHERE";
            name = :x, module_id = :Central, prior = LogNormal(0.0, 1.0))
    end

    @testset "Governing choices are enumerable from the loaded parameters" begin
        pgk3 = load_parameter(both, "kcat_fwd_PGK3";
                              name = :kcat_PGK3, module_id = :Nucleotide,
                              prior = LogNormal(0.0, 1.0),
                              governing = "nucleotide_balanced")
        eno = load_parameter(both, "kcat_fwd_ENO";
                             name = :kcat_ENO, module_id = :Central,
                             prior = LogNormal(0.0, 1.0))

        choices = governing_choices([pgk3, eno])
        @test length(choices) == 1              # only PGK3 had an alternative
        name, chosen, rejected = choices[1]
        @test name == :kcat_PGK3
        @test chosen == "nucleotide_balanced"
        @test rejected == ["central_balanced"]
    end

    @testset "Asserted priors are distinguishable from inherited ones" begin
        pts = load_parameter(both, "k_pts_GLCpts0_fwd";
                             name = :k_GLCpts0_fwd, module_id = :Transport,
                             prior = LogNormal(0.0, 2.0))
        eno = load_parameter(both, "kcat_fwd_ENO";
                             name = :kcat_ENO, module_id = :Central,
                             prior = LogNormal(0.0, 1.0))
        # G3P sits at the prior median at prior width in the central file — and
        # at exactly that value in the registry, so the agreement check passes.
        g3p = load_parameter(both, "conc_M_g3p_c";
                             name = :g3p0, module_id = :Central,
                             prior = LogNormal(0.0, 1.0),
                             role = :initial_condition)

        params = [pts, eno, g3p]

        # The PTS mass-action constants carry no quantified uncertainty in the
        # source, so any prior on them is this project's.
        @test [p.name for p in asserted_prior_params(params)] == [:k_GLCpts0_fwd]
        # The balanced glycolytic parameter is not asserted.
        @test informedness(eno) == :balanced
        # And the prior-median row is flagged as uninformed.
        @test [p.name for p in uninformed_params(params)] == [:g3p0]
    end

    @testset "The registry and the loader cannot disagree on a concentration" begin
        # GTP is a nucleotide-module species: the registry records 1.6627 from
        # the nucleotide file. Importing the central file's 0.1 prior default —
        # the cross-file error of record — now fails instead of passing with
        # clean provenance.
        err = caught() do
            load_parameter(both, "conc_M_gtp_c";
                           name = :gtp0, module_id = :Nucleotide,
                           prior = LogNormal(0.0, 1.0),
                           role = :initial_condition,
                           governing = "central_balanced")
        end
        @test err isa ArgumentError
        @test occursin("0.1", err.msg)
        @test occursin("1.6627", err.msg)
        @test occursin("M_gtp_c", err.msg)

        # The governing file the registry agrees with loads cleanly.
        gtp = load_parameter(both, "conc_M_gtp_c";
                             name = :gtp0, module_id = :Nucleotide,
                             prior = LogNormal(0.0, 1.0),
                             role = :initial_condition,
                             governing = "nucleotide_balanced")
        @test gtp.value == 1.6627
        @test informedness(gtp) == :balanced

        # Identifiers that name no registry species are unconstrained.
        @test load_parameter(both, "kcat_fwd_ENO";
                             name = :kcat_ENO, module_id = :Central,
                             prior = LogNormal(0.0, 1.0)).value == 62.18
    end

    @testset "A governing declaration binds even with a single holder" begin
        # Only the central table holds ENO. Declaring the nucleotide file as
        # governing means the table set and the declaration disagree — a stale
        # table, not a resolvable ambiguity.
        err = caught() do
            load_parameter(both, "kcat_fwd_ENO";
                           name = :kcat_ENO, module_id = :Central,
                           prior = LogNormal(0.0, 1.0),
                           governing = "nucleotide_balanced")
        end
        @test err isa ArgumentError
        @test occursin("kcat_fwd_ENO", err.msg)
        @test occursin("nucleotide_balanced", err.msg)
        @test occursin("central_balanced", err.msg)

        # The matching declaration is redundant but consistent, so it loads.
        p = load_parameter(both, "kcat_fwd_ENO";
                           name = :kcat_ENO, module_id = :Central,
                           prior = LogNormal(0.0, 1.0),
                           governing = "central_balanced")
        @test p.value == 62.18
    end

    @testset "A truncated row is corruption, not something to skip" begin
        # Silently dropping a short row could collapse a two-file ambiguity to
        # one holder and bypass the governs machinery entirely.
        mktempdir() do dir
            path = joinpath(dir, "truncated.tsv")
            write(path, "!ID\t!Mode\t!GeometricStd\nkcat_fwd_PGK3\n")
            err = caught(() -> read_source_table(path; file = "truncated"))
            @test err isa ArgumentError
            @test occursin("truncated", err.msg)
        end
    end

    @testset "A row truncated before the uncertainty columns is also corruption" begin
        # Cut short before Informedness, the row would otherwise silently
        # re-derive an informedness it declared — the declaration rot the
        # declared-informedness check exists to reject.
        mktempdir() do dir
            path = joinpath(dir, "short.tsv")
            write(path, "!ID\t!Mode\t!GeometricStd\t!Informedness\nkcat_fwd_X\t1.0\t20.0\n")
            @test_throws ArgumentError read_source_table(path; file = "short")

            # And cut short before GeometricStd, it would silently classify
            # :asserted.
            write(path, "!ID\t!Mode\t!GeometricStd\nkcat_fwd_X\t1.0\n")
            @test_throws ArgumentError read_source_table(path; file = "short")
        end
    end

    @testset "A non-empty unparseable uncertainty is rejected, not misfiled" begin
        mktempdir() do dir
            path = joinpath(dir, "gstd.tsv")
            # "N/A" would otherwise classify :asserted — corruption misfiled as
            # a point value.
            write(path, "!ID\t!Mode\t!GeometricStd\nkcat_fwd_X\t1.0\tN/A\n")
            err = caught(() -> read_source_table(path; file = "gstd"))
            @test err isa ArgumentError
            @test occursin("N/A", err.msg)

            # NaN compares false against the prior width, which would report an
            # unquantified value as informed.
            write(path, "!ID\t!Mode\t!GeometricStd\nkcat_fwd_X\t1.0\tNaN\n")
            @test_throws ArgumentError read_source_table(path; file = "gstd")

            # An empty cell legitimately means the source asserts a point value.
            write(path, "!ID\t!Mode\t!GeometricStd\nkcat_fwd_X\t1.0\t\n")
            t = read_source_table(path; file = "gstd")
            @test t.informedness["kcat_fwd_X"] == :asserted
        end
    end

    @testset "The default logical name is the extension-free basename" begin
        # What provenance records and governing declarations match is the
        # registry's logical source name, not a filename with an extension.
        t = read_source_table(joinpath(FIXTURES, "central_balanced.tsv"))
        @test t.file == "central_balanced"
    end

    @testset "Registry agreement tolerates the 4-decimal transcription" begin
        # The registry transcribes 1.6627; a real balanced table carries full
        # precision. Agreement means equal up to the transcription, not
        # bit-for-bit.
        mktempdir() do dir
            path = joinpath(dir, "fullprec.tsv")
            write(path, "!ID\t!Mode\t!GeometricStd\nconc_M_gtp_c\t1.66271\t1.32\n")
            t = read_source_table(path; file = "nucleotide_balanced")
            p = load_parameter([t], "conc_M_gtp_c";
                               name = :gtp0, module_id = :Nucleotide,
                               prior = LogNormal(0.0, 1.0),
                               role = :initial_condition)
            @test p.value == 1.66271

            # A genuinely different number is still two copies of one number.
            write(path, "!ID\t!Mode\t!GeometricStd\nconc_M_gtp_c\t1.67\t1.32\n")
            stale = read_source_table(path; file = "nucleotide_balanced")
            @test_throws ArgumentError load_parameter([stale], "conc_M_gtp_c";
                name = :gtp0, module_id = :Nucleotide,
                prior = LogNormal(0.0, 1.0), role = :initial_condition)
        end
    end

    @testset "A misspelled declared informedness is rejected, not re-derived" begin
        mktempdir() do dir
            path = joinpath(dir, "declared.tsv")
            write(path, "!ID\t!Mode\t!Informedness\nkcat_fwd_X\t1.0\tblanced\n")
            err = caught(() -> read_source_table(path; file = "declared"))
            @test err isa ArgumentError
            @test occursin("blanced", err.msg)

            # A valid declaration is honoured, and an empty cell falls back to
            # derivation.
            write(path, "!ID\t!Mode\t!Informedness\nkcat_fwd_X\t1.0\tasserted\nkcat_fwd_Y\t1.0\t\n")
            t = read_source_table(path; file = "declared")
            @test t.informedness["kcat_fwd_X"] == :asserted
            @test t.informedness["kcat_fwd_Y"] == :asserted   # no gstd column
        end
    end

    @testset "A directly constructed parameter reports provenance as absent" begin
        p = InferParameter(1.0, Normal(0, 1), false, :k, :mod, :rate)
        @test provenance_of(p) === nothing
        @test source_file(p) === nothing
        @test informedness(p) == :not_imported
    end

    @testset "Provenance conflicts are found across modules" begin
        from_central = load_parameter(both, "kcat_fwd_PGK3";
                                      name = :kcat_PGK3, module_id = :Central,
                                      prior = LogNormal(0.0, 1.0),
                                      governing = "central_balanced")
        from_nucleotide = load_parameter(both, "kcat_fwd_PGK3";
                                         name = :kcat_PGK3, module_id = :Nucleotide,
                                         prior = LogNormal(0.0, 1.0),
                                         governing = "nucleotide_balanced")
        conflicts = provenance_conflicts([from_central, from_nucleotide])
        @test length(conflicts) == 1
        @test conflicts[1][1] == :kcat_PGK3
        @test conflicts[1][2] == ["central_balanced", "nucleotide_balanced"]

        # One source only: no conflict.
        @test isempty(provenance_conflicts([from_central]))
    end

    @testset "Provenance survives deduplication" begin
        tagged = load_parameter(both, "kcat_fwd_ENO";
                                name = :kcat_ENO, module_id = :Central,
                                prior = LogNormal(0.0, 1.0))
        plain = InferParameter(1.0, Normal(0, 1), false, :k_other, :Other, :rate)

        deduped = unique_params([tagged, plain, tagged])
        @test length(deduped) == 2
        survivor = only(filter(p -> p.name == :kcat_ENO, deduped))
        @test source_file(survivor) == "central_balanced"
        @test provenance_of(survivor).identifier == "kcat_fwd_ENO"

        # A provenance-less copy seen first does not shadow a tagged one seen
        # later: which copy survives must not depend on composition order.
        untracked = InferParameter(62.18, LogNormal(0.0, 1.0), false,
                                   :kcat_ENO, :Other, :rate)
        for ordering in ([tagged, untracked], [untracked, tagged])
            kept = only(unique_params(ordering))
            @test source_file(kept) == "central_balanced"
        end

        # And the role filters carry it through too.
        @test source_file(only(rate_params([tagged]))) == "central_balanced"
        @test source_file(only(free_params([tagged]))) == "central_balanced"
    end

    @testset "Positional ParameterSource construction cannot bypass validation" begin
        # Validation lives in the inner constructor, as for the edge kinds: an
        # off-vocabulary informedness would silently drop the parameter from
        # every provenance report.
        @test_throws ArgumentError ParameterSource(
            "central_balanced", nothing, nothing, :ballanced,
            Pair{String, Float64}[])
    end

    @testset "An initial condition outside the conc_ convention warns" begin
        # The registry-agreement check keys on the conc_<species> prefix, so an
        # initial condition named otherwise forfeits the check silently — the
        # role signal makes the forfeiture visible.
        mktempdir() do dir
            path = joinpath(dir, "renamed.tsv")
            write(path, "!ID\t!Mode\t!GeometricStd\nIC_M_gtp_c\t1.6627\t1.32\n")
            t = read_source_table(path; file = "nucleotide_balanced")
            p = @test_logs (:warn, r"conc_") load_parameter(
                [t], "IC_M_gtp_c";
                name = :gtp0, module_id = :Nucleotide,
                prior = LogNormal(0.0, 1.0), role = :initial_condition)
            @test p.value == 1.6627

            # A rate import outside the convention is not an initial condition
            # and warns about nothing.
            write(path, "!ID\t!Mode\t!GeometricStd\nkcat_fwd_X\t1.0\t1.3\n")
            t2 = read_source_table(path; file = "misc")
            @test_logs min_level=Base.CoreLogging.Warn load_parameter(
                [t2], "kcat_fwd_X";
                name = :k_x, module_id = :Central, prior = LogNormal(0.0, 1.0))
        end
    end
end
