struct ObservedData
    times::Vector{Float64}
    observations::Matrix{Float64}  # species x timepoints
    species::Vector{Symbol}
end

struct PosteriorPredictive
    solutions::Vector
    times::Vector{Float64}
    species::Vector{Symbol}
end
