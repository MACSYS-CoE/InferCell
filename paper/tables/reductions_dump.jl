# Prints the declared reductions of the reduced syn3A model and of M0, from
# `reduction_declarations`, as tab-separated composition, category, subject, description.
# make_modelform.py turns the output, reductions_raw.tsv, into the paper's table's counts.
#
# Run from the repository root:
#   julia --project paper/tables/reductions_dump.jl > paper/tables/reductions_raw.tsv
using InferCell

println("composition\tcategory\tsubject\tdescription")
for (name, models, driver) in (("reduced", corea_models(), build_corea()),
                               ("m0", m0_models(), build_m0()))
    for l in reduction_declarations(models, driver)
        println(join((name, l.category, l.subject, replace(l.description, r"\s+" => " ")), '\t'))
    end
end
