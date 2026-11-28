# KIM-56: the named kinetic constants that the override hook exposes, as a table.
# Writes OUT/param_names.csv (index,name,module,reaction,key,default) and prints the
# counts the result note quotes. Exits non-zero if a count is not as expected.
#
# Environment: OUT (output directory).

using SpeedyCell.Cell
using SpeedyCell.ODE: K_NAMES
using Printf

function main()
    out = get(ENV, "OUT", "")
    isempty(out) && error("set OUT")
    mkpath(out)
    names = param_names()
    defaults = param_defaults()
    ode = Set(K_NAMES)
    rows = map(enumerate(names)) do (i, n)
        if n in ode
            j = findlast('.', n)
            (i, n, "ODE", n[1:j-1], n[j+1:end], defaults[n])
        else
            (i, n, "CME", "", "", defaults[n])
        end
    end
    open(joinpath(out, "param_names.csv"), "w") do io
        println(io, "index,name,module,reaction,key,default")
        for (i, n, mod, r, k, d) in rows
            q = occursin(',', n) || occursin('"', n) ? "\"" * replace(n, "\"" => "\"\"") * "\"" : n
            println(io, i, ',', q, ',', mod, ',', r, ',', k, ',', repr(d))
        end
    end
    nODE = count(r -> r[3] == "ODE", rows)
    nCME = count(r -> r[3] == "CME", rows)
    vals = [r[6] for r in rows]
    nzero = count(==(0.0), vals)
    nbad = count(x -> !isfinite(x) || x < 0, vals)
    keys = Dict{String,Int}()
    for r in rows
        r[3] == "ODE" && (keys[r[5]] = get(keys, r[5], 0) + 1)
    end
    println("named constants: ", length(names), " (unique ", length(unique(names)), "), ODE ", nODE,
            ", CME ", nCME, "; ODE reactions with a named constant: ",
            length(unique(r[4] for r in rows if r[3] == "ODE")))
    println("defaults: zero ", nzero, ", negative or non-finite ", nbad,
            @sprintf(", range %.3g to %.3g", minimum(vals), maximum(vals)))
    println("ODE keys: ", join(["$k $v" for (k, v) in sort(collect(keys); by=x -> (-x[2], x[1]))], ", "))
    println("ODE-first order: ", all(r -> r[3] == "ODE", rows[1:nODE]))
    ok = length(names) == 838 && length(unique(names)) == 838 && nODE == 805 && nCME == 33 &&
         nzero == 0 && nbad == 0 && all(r -> r[3] == "ODE", rows[1:nODE])
    println(ok ? "TABLE PASS" : "TABLE FAIL")
    ok || exit(1)
end

main()
