abstract type StemCell end
abstract type StemCellModel end

mutable struct StemCell <: StemCell
    r::Float64
    λ::Float64
    Δ::Float64
    μd::Float64
    μp::Float64
    mutationsd::Array{Int64, 1}
    mutationsp::Array{Int64, 1}
end

mutable struct StemCellModel <: StemCellModel
    Ns::Int64
    tend::Float64
    μd::Float64
    μp::Float64
    r::Float64
    λ::Float64
    Δ::Float64
    Δmut::Float64
    ρ::Float64
    s::Float64

    StemCellModel(Ns::Int64;
    tend = 10.0,
    μd = 0.01,
    μp = 1.0,
    r = 0.5,
    λ = 1.0,
    Δ = 0.0,
    Δmut = 0.0,
    ρ = 1.0,
    s = 0.0) = new(Ns, tend, μd, μp, r, λ, Δ, Δmut, ρ, s)
end

struct StemCellResults
    stemcells::Array{StemCell, 1}
    SM::StemCellModel
    clonesize::DataFrame
    mutationfrequencies_d::Array{Float64, 1}
    mutationfrequencies_p::Array{Float64, 1}
    mutationsize_d::Array{Float64, 1}
    mutationsize_p::Array{Float64, 1}
    hitchikers::Array{Int64, 1}

    StemCellResults(stemcells, SM) =
    new(stemcells, SM, clonesize(stemcells, SM),
    mutationfrequencies(stemcells, SM)[1],
    mutationfrequencies(stemcells, SM)[2],
    mutationsize(stemcells, SM)[1],
    mutationsize(stemcells, SM)[2],
    mutationsize(stemcells, SM)[3])

end

function copystemcell(sc::StemCell)
    StemCell(copy(sc.r), copy(sc.λ) , copy(sc.Δ), copy(sc.μd), copy(sc.μp), copy(sc.mutationsd), copy(sc.mutationsp))
end


function createstemcellpool(SM::StemCellModel)
    stemcellpool = StemCell[]
    for _ in 1:SM.Ns
        push!(stemcellpool,
        StemCell(SM.r, SM.λ, SM.Δ, SM.μd, SM.μp, Int64[], Int64[]))
    end
    return stemcellpool
end

function newmutations(sc::StemCell, mutID, SM::StemCellModel; onedriver = true)

    #function to add new mutations to cells based on μ
    #numbermutations = μ
    #println("ondriver is $onedriver")
    numbermutationsp = rand(Poisson(sc.μp))
    numbermutationsd = rand(Poisson(sc.μd))
    if length(sc.mutationsd) > 0 && onedriver
        numbermutationsd = 0
        numbermutationsp = 0
    end
    append!(sc.mutationsp, mutID:(mutID + numbermutationsp - 1))
    mutID = mutID + numbermutationsp
    append!(sc.mutationsd, mutID:(mutID + numbermutationsd - 1))
    mutID = mutID + numbermutationsd

    if numbermutationsd > 0
        sc.Δ = SM.Δmut
        sc.λ = (1 + SM.s) * sc.λ
    end

    return sc, mutID
end

function cellturnover(scpool, randomsc, mutID, SM; onedriver = true)
    r = rand()
    if r < ((1 + scpool[randomsc].Δ) / 2)
        push!(scpool, copystemcell(scpool[randomsc]))
        scpool[end], mutID =
            newmutations(scpool[end], mutID, SM, onedriver = onedriver)
    else
        deleteat!(scpool, randomsc)
    end
    return scpool, mutID
end

function cellturnover(scpool, randomsc; onedriver = true)
    r = rand()
    if r < ((1 + scpool[randomsc].Δ) / 2)
        push!(scpool, copystemcell(scpool[randomsc]))
        rdaughter = rand([randomsc, length(scells)])
        scpool[rdaughter], mutID =
            newmutations(scpool[rdaughter], mutID, SM, onedriver = onedriver)
    else
        deleteat!(scpool, randomsc)
    end
    return scpool, mutID
end

function runsimulation(SM::StemCellModel; progress = false, restart = false, finish = "time", onedriver = true)
    scpool = createstemcellpool(SM)
    t = 0.0
    Rmax = SM.λ * (SM.s + 1) * SM.r * 2
    N = length(scpool)
    mutID = 1
    Nmax = round(N*exp(Rmax * SM.tend))
    Ndivs = 0

    if progress == true
        p = ProgressThresh(0.0, "Time left to simulate:")
    end

    if finish == "time"
        while t < SM.tend
            randomsc = rand(1:length(scpool))
            u1 = rand(Uniform(0, Rmax))
            if u1 < (scpool[randomsc].λ * scpool[randomsc].r * 2)
                #println("divide")
                scpool, mutID = cellturnover(scpool, randomsc, mutID, SM, onedriver = onedriver)
                Δt = (1 / (Rmax * N)) * log(1/rand())
                t += Δt
                N = length(scpool)
                Ndivs += 1
            else
                #println("die")
                Δt = (1 / (Rmax * N)) * log(1/rand())
                t += Δt
                N = length(scpool)
            end
            if progress == true
                ProgressMeter.update!(p, SM.tend - t)
            end
            if N == 0
                #@warn "Population died out"
                if restart == true
                    scpool = createstemcellpool(SM)
                    t = 0.0
                    Rmax = SM.λ * (SM.s + 1) * SM.r * 2
                    N = length(scpool)
                    mutID = 1
                    Nmax = round(N*exp(Rmax * SM.tend))
                else
                    return StemCellResults(scpool, SM)
                end
            end
        end
    else
        #println("Nmax is $(Nmax)")
        while N < Nmax
            randomsc = rand(1:length(scpool))
            u1 = rand(Uniform(0, Rmax))
            if u1 < (scpool[randomsc].λ * scpool[randomsc].r * 2)
                #println("divide")
                scpool, mutID = cellturnover(scpool, randomsc, mutID, SM, onedriver = onedriver)
                Δt = (1 / (Rmax * N)) * log(1/rand())
                t += Δt
                N = length(scpool)
                Ndivs += 1
            else
                #println("die")
                Δt = (1 / (Rmax * N)) * log(1/rand())
                t += Δt
                N = length(scpool)
            end
            if progress == true
                ProgressMeter.update!(p, SM.tend - t)
            end
            if N == 0
                #@warn "Population died out"
                if restart == true
                    scpool = createstemcellpool(SM)
                    t = 0.0
                    Rmax = SM.λ * SM.r * 2
                    N = length(scpool)
                    mutID = 1
                    Nmax = round(N*exp(Rmax * SM.tend))
                else
                    return StemCellResults(scpool, SM)
                end
            end
        end
    end
    return StemCellResults(scpool, SM)
end

function cellsconvert(scells)
    #convert from array of cell types to one array with mutations and one array with cell fitness
    mutationsp = Int64[]
    mutationsd = Int64[]
    hitchiker = Int64[]
    for i in 1:length(scells)
        driver = Int64(length(scells[i].mutationsd) > 0)
        append!(mutationsp, scells[i].mutationsp)
        append!(mutationsd, scells[i].mutationsd)
        append!(hitchiker, fill(driver, length(scells[i].mutationsp)))
    end
    return mutationsp, mutationsd, hitchiker
end

Pn(r, λ, t, n) = (1 / log.(r * λ * t)) .* exp.(-n / (r * λ * t)) .* (1 ./ n)
function Pn(r, λ, t, n, Δ)
    Ncp = ((1 + Δ) * exp(2 * r * λ * Δ * t) - (1 - Δ)) / (2 * Δ)
    return (1 ./ n) .* (1 / log(Ncp)) .* exp.(-n ./ Ncp)
end

function mutationfrequencies(scells, SM)
    mp, md, h = cellsconvert(scells)
    m = [mp; md]
    if length(mp) == 0
        mutationfrequency_passengers = []
        hcounts = []
    else
        mutationfrequency_passengers = StatsBase.counts(mp, 1:maximum(mp)) ./ length(scells)
        mcounts = StatsBase.counts(mp, 1:maximum(mp))
        hcounts = fill(0, length(mcounts))
        hcounts[mp[findall(x -> x > 0, h)]] .= 1
        idxkeep = findall(x -> x > 0.0, mutationfrequency_passengers)
        mutationfrequency_passengers = mutationfrequency_passengers[idxkeep]
        hcounts = hcounts[idxkeep]
    end

    if length(md) == 0
        mutationfrequency_drivers = []
    else
        mutationfrequency_drivers = StatsBase.counts(md, 1:maximum(md)) ./ length(scells)
        filter!(x -> x > 0.0, mutationfrequency_drivers)
    end
    return mutationfrequency_drivers, mutationfrequency_passengers, hcounts
end

function mutationsize(scells, SM)
    mp, md, h = cellsconvert(scells)
    m = [mp; md]
    if length(mp) == 0
        mutationfrequency_passengers = []
        hcounts = []
    else
        mutationfrequency_passengers = StatsBase.counts(mp, 1:maximum(mp))
        mcounts = StatsBase.counts(mp, 1:maximum(mp))
        hcounts = fill(0, length(mcounts))
        hcounts[mp[findall(x -> x >0, h)]] .= 1
        idxkeep = findall(x -> x > 0.0, mutationfrequency_passengers)
        mutationfrequency_passengers = mutationfrequency_passengers[idxkeep]
        hcounts = hcounts[idxkeep]
    end

    if length(md) == 0
        mutationfrequency_drivers = []
    else
        mutationfrequency_drivers = StatsBase.counts(md, 1:maximum(md))
        filter!(x -> x > 0.0, mutationfrequency_drivers)
    end
    return mutationfrequency_drivers, mutationfrequency_passengers, hcounts
end


function clonesize(scells, SM)
    mp, md, h = cellsconvert(scells)
    m = [mp; md]
    if length(m) == 0
        return DataFrame()
    end
    clonefrequency = StatsBase.counts(sort(StatsBase.counts(m, 1:maximum(m))), 1:length(scells))
    cloneprobability = clonefrequency / sum(clonefrequency)
    DF = DataFrame(
    n = 1:length(scells),
    A = 1:length(scells) ./ SM.ρ,
    Cn = clonefrequency,
    Pn = cloneprobability
    )
    DF[:Pntheory] = Pn(SM.r, SM.λ, SM.tend, DF[:n])
    DF[:Pntheoryselection] = Pn(SM.r, SM.λ, SM.tend, DF[:n], SM.Δ)
    return DF
end

function averageclonesize(results::StemCellResults)
    avesize = sum(results.clonesize[:Pn] .* results.clonesize[:n])
    Ncp = (results.SM.r * results.SM.λ * results.SM.tend)
    neuttheory = (Ncp - 1) / (log(Ncp))
    Ncp = ((1 + results.SM.Δmut) * exp(2 * results.SM.r * results.SM.λ * results.SM.Δmut * results.SM.tend) - (1 - results.SM.Δmut)) / (2 * results.SM.Δmut)
    seltheory = (Ncp - 1) / (log(Ncp))
    return DataFrame(data = avesize, theoryneut = neuttheory,
    seltheory = seltheory)
end

function averageclonesize2(results::StemCellResults)
    avesize = sum(results.clonesize[:Pn] .* results.clonesize[:n])
    Ncp = (results.SM.r * results.SM.λ * results.SM.tend)

    neuttheory = 1 / (log(Ncp) * (exp(1/Ncp) - 1))
    Ncp = ((1 + results.SM.Δ) * exp(2 * results.SM.r * results.SM.Δ * results.SM.λ * results.SM.tend) - (1 - results.SM.Δ)) / (2 * results.SM.Δ)
    seltheory = 1 / (log(Ncp) * (exp(1/Ncp) - 1))
    return DataFrame(data = avesize, theoryneut = neuttheory,
    seltheory = seltheory)
end

function show(io::IO, results::StemCellResults)
    println("######################################")
    println("Input parameters:")
    println("   Number of Stem Cells: $(results.SM.Ns)")
    println("   rλ: $(results.SM.r * results.SM.λ)")
    println("   Δ: $(results.SM.Δ)")
    println("   Δmut: $(results.SM.Δmut)")
    println("   μd: $(results.SM.μd)")
    println("   μp: $(results.SM.μp)")
    println("   ρ: $(results.SM.ρ)")
    println("   t: $(results.SM.tend)")
    println()
    println("######################################")
    println("Results:")
    println("   Number of cells: $(length(results.stemcells))")
    println("   % cells with driver: $(sum(map(x -> length(x.mutationsd), results.stemcells) .> 0) / length(results.stemcells))")
    println("   Mean bias: $(mean(map(x -> x.Δ, results.stemcells)))")
    println("   Mean driver clone size: $(mean(results.mutationsize_d))")
    if length(results.mutationsize_d) > 5
        println(UnicodePlots.histogram(results.mutationsize_d, nbins = 20))
        println()
    end
    println("   Mean passenger clone size: $(mean(results.mutationsize_p))")
    if length(results.mutationsize_p) > 5
        println(UnicodePlots.histogram(results.mutationsize_p, nbins = 20))
        println()
    end
    println("   Mean driver clone frequency: $(mean(results.mutationfrequencies_d))")
    println("   Mean passenger clone frequency: $(mean(results.mutationfrequencies_p))")
    println()
    println("######################################")
end
