using Test
using InferCell
using Random
using LinearAlgebra

# Spec §11 phase 14b: the derivative, ensemble and census checks on the
# assembled Core A′, each with a mutation showing it can fail.
#
# The full-scale runs are in `dev/scripts/corea_validation_14b.jl` (check 1b's
# ensemble over a cycle, check 6, task 9.9's rerun and F5) and
# `dev/scripts/corea_census.jl` (check 7 over seeds and 200 prior draws).
# Their result files are the numbers of record. This file asserts each check's
# mechanism and its mutation at a size the suite can afford.

const B_SEED = 1410
const B_SHORT = 600

@testset "Phase 14b: the derivative, ensemble and census checks" begin

    # ------------------------------------------------------------------
    @testset "14.8 check 8: the summation theorems on the frozen variant" begin
        cq = frozen_control_problem()
        @test length(cq.multipliers) == 24          # 22 ODE reactions and two demands

        # The conserved combinations: every moiety of 14a that the frozen
        # variant keeps lies in the left null space of its stoichiometry.
        # Carbon is not among them: glucose enters and lactate leaves.
        cm = conservation_matrices(cq)
        names = cq.names[cq.internal]
        inL(pairs) = (v = zeros(length(names));
                      for (s, w) in pairs; v[findfirst(==(s), names)] = w; end;
                      norm(v - cm.L' * (cm.L * v)) / norm(v))
        for m in corea_moieties(carbon = false)
            @test inL(m.ode) < 1e-12
        end
        @test inL((:M_trna_c => 1, :M_trna_chg_c => 1)) < 1e-12
        # And one more, which no 14a check asserts (see the result file).
        @test size(cm.L, 1) == 10

        ss = steady_state(cq)
        @test ss.residual <= 1e-12 * ss.scale
        @test all(<(0), real.(ss.eigenvalues))       # stable in its conserved class
        @test all(>=(0), ss.x)

        cc = control_coefficients(cq, ss)
        # GK1's GMP has no source once expression is frozen, so it carries no
        # steady flux, and its flux has no control coefficient.
        @test cc.zero_fluxes == [:R_GK1]
        dev = assert_summation(cc)
        @test dev.ccc <= 1e-6
        @test dev.fcc <= 1e-6

        # The mutation spec §3 names: without k_chg's multiplier the identities
        # fail, and the check names the quantity.
        mut = omit_multiplier(cc, :R_charging)
        err = try
            assert_summation(mut); nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Check 8 fails", err.msg)
        live = [k for k in eachindex(mut.multipliers) if !(mut.multipliers[k] in mut.zero_fluxes)]
        @test maximum(abs.(mut.fcc_sum[live] .- 1)) > 1e-2
        # The charging flux's own row loses exactly its elasticity to itself.
        k = findfirst(==(:R_charging), cc.multipliers)
        @test cc.fcc[k, k] > 0.01

        # Gene-level coefficients sum a gene's reactions; they satisfy the
        # concentration identity only with the non-gene multipliers added back.
        Cg, genes = group_coefficients(cc.ccc, cc.multipliers, frozen_gene_groups())
        @test length(genes) == 13
        rest = [j for (j, m) in enumerate(cc.multipliers)
                if !any(m in last(g) for g in frozen_gene_groups())]
        @test maximum(abs.(vec(sum(Cg; dims = 2)) .+ vec(sum(cc.ccc[:, rest]; dims = 2)))) <= 1e-6

        # The problem is not mutated by the analysis.
        @test cq.p == frozen_control_problem().p
    end

end
