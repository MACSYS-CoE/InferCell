# M0's after-run audit (spec §12, KIM-40): every record's ran_on, commit and SMOKE.
using Serialization, Dates, Random
using InferCell

const OUT = ARGS[1]
const COMMIT = "ed179c28c561a54ef931e94a54c4e8b6e9cc99d2"
const TURIN = ("dave148", "dave149", "gina18")

function audit()
    files = String[joinpath(OUT, "manifest.jls"), joinpath(OUT, "score.jls")]
    for d in ("work", "work_plain", "finish")
        append!(files, filter(endswith(".jls"), readdir(joinpath(OUT, d); join = true)))
    end
    cpus = Dict{String,Int}(); hosts = Dict{String,Int}(); jobs = Dict{String,Int}()
    bad = String[]
    for f in files
        try
            x = deserialize(f)
            r = x.ran_on
            cpus[r.cpu_name] = get(cpus, r.cpu_name, 0) + 1
            hosts[r.host] = get(hosts, r.host, 0) + 1
            jobs[r.array_job] = get(jobs, r.array_job, 0) + 1
            r.cpu_name == "znver3" || push!(bad, "$f: cpu_name $(r.cpu_name)")
            r.host in TURIN && push!(bad, "$f: host $(r.host)")
            x.commit_full == COMMIT || push!(bad, "$f: commit $(x.commit_full)")
            x.smoke == false || push!(bad, "$f: smoke $(x.smoke)")
            if startswith(basename(dirname(f)), "work")
                i = parse(Int, basename(f)[1:4])
                r.array_task == string(i) || push!(bad, "$f: array_task $(r.array_task)")
            end
        catch e
            push!(bad, "$f: unreadable: $(sprint(showerror, e))")
        end
    end
    println("records read: ", length(files))
    println("cpu_name: ", cpus)
    println("array_job: ", jobs)
    println("nodes: ", length(hosts), "; gina: ", count(startswith("gina"), keys(hosts)),
            "; on a turin node: ", count(in(TURIN), keys(hosts)))
    println("problems: ", length(bad))
    foreach(println, first(bad, 50))
    return isempty(bad)
end

audit() || exit(1)
