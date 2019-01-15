mutable struct SkinStemCell2D <: StemCell
    r::Float64
    λ::Float64
    Δ::Float64
    μd::Float64
    μp::Float64
    mutationsd::Array{Int64, 1}
    mutationsp::Array{Int64, 1}
    x::Int64
    y::Int64
end

mutable struct SkinStemCellModel2D <: StemCellModel
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
    Nx::Int64
    Ny::Int64

    SkinStemCellModel2D(Ns::Int64;
    tend = 10.0,
    μd = 0.01,
    μp = 1.0,
    r = 0.5,
    λ = 1.0,
    Δ = 0.0,
    Δmut = 0.0,
    ρ = 1.0,
    s = 0.0,
    Nx = 100,
    Ny = 100) = new(Ns, tend, μd, μp, r, λ, Δ, Δmut, ρ, s, Nx, Ny)
end

struct SkinStemCellResults2D
    stemcells::Array{SkinStemCell2D, 1}
    SM::SkinStemCellModel2D

    SkinStemCellResults2D(stemcells, SM) =
    new(stemcells, SM)
end

function copyskinstemcell(sc::SkinStemCell2D)
    SkinStemCell2D(copy(sc.r), copy(sc.λ) , copy(sc.Δ), copy(sc.μd), copy(sc.μp), copy(sc.mutationsd), copy(sc.mutationsp),
    copy(sc.x), copy(sc.y))
end

function createstemcellpool(SM::SkinStemCellModel2D)
    midpoint = [round(SM.Nx ./ 2), round(SM.Nx ./ 2)]
    radius = 1.5 * sqrt(SM.Ns / π)
    stemcellpool = SkinStemCell2D[]
    xycoords = String[]
    scids = DefaultDict(0, Dict((midpoint[1], midpoint[2]) => 1))
    push!(stemcellpool,
    SkinStemCell2D(SM.r, SM.λ, SM.Δ, SM.μd, SM.μp, Int64[], Int64[], midpoint[1], midpoint[2]))
    Nstemcells = 1
    while Nstemcells < SM.Ns
        r = radius * sqrt(rand())
        θ = rand() * 2 * π
        x = round(midpoint[1] + r * cos(θ))
        y = round(midpoint[2] + r * sin(θ))
        if scids[(x,y)] == 0
            push!(stemcellpool,
            SkinStemCell2D(SM.r, SM.λ, SM.Δ, SM.μd, SM.μp, Int64[], Int64[], x, y))
            Nstemcells += 1
            push!(scids, (x, y) => Nstemcells)
        end
    end
    return stemcellpool, scids
end

function placestemcell(scpool, randomsc, scids)
    direction = rand([-1, 1], 2)
    newx = scpool[randomsc].x + direction[1]
    newy = scpool[randomsc].y + direction[2]
    println("Current stem cell: ($(scpool[randomsc].x), $(scpool[randomsc].y), index:$(randomsc))")
    println("New stem cell: ($(newx), $(newy))")
    println("Size of scpool: $(length(scpool)),
        id of new sc: $(scids[(newx,newy)])")
    if scids[(newx,newy)] == 0
        newstemcell = copyskinstemcell(scpool[randomsc])
        newstemcell.x = newx
        newstemcell.y = newy
        push!(scpool, newstemcell)
        push!(scids, (newx, newy) => length(scpool))
        return scpool, scids
    else
        deleteat!(scpool, scids[(newx,newy)])
        newstemcell = copyskinstemcell(scpool[randomsc])
        newstemcell.x = newx
        newstemcell.y = newy
        push!(scpool, newstemcell)
        push!(scids, (newx, newy) => length(scpool))
        return scpool, scids
    end
end

function cellturnover(scpool, randomsc, mutID, SM, scids; onedriver = true)
    r = rand()
    if r < ((1 + scpool[randomsc].Δ) / 2)
        scpool, scids = placestemcell(scpool, randomsc, scids)
        scpool[end], mutID =
            newmutations(scpool[end], mutID, SM, onedriver = onedriver)
    else
        deleteat!(scpool, randomsc)
    end
    return scpool, mutID, scids
end

function runsimulation(SM::SkinStemCellModel2D; progress = false, restart = false, finish = "time", onedriver = true)
    scpool, scids = createstemcellpool(SM)
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
                scpool, mutID, scids = cellturnover(scpool, randomsc,
                                        mutID, SM, scids,
                                        onedriver = onedriver)
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
                update!(p, SM.tend - t)
            end
            if N == 0
                #@warn "Population died out"
                if restart == true
                    scpool, scids = createstemcellpool(SM)
                    t = 0.0
                    Rmax = SM.λ * (SM.s + 1) * SM.r * 2
                    N = length(scpool)
                    mutID = 1
                    Nmax = round(N*exp(Rmax * SM.tend))
                else
                    return SkinStemCellResults2D(scpool, SM)
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
                scpool, mutID, scids = cellturnover(scpool, randomsc,
                                        mutID, SM, scids,
                                        onedriver = onedriver)
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
                update!(p, SM.tend - t)
            end
            if N == 0
                #@warn "Population died out"
                if restart == true
                    scpool, scids = createstemcellpool(SM)
                    t = 0.0
                    Rmax = SM.λ * SM.r * 2
                    N = length(scpool)
                    mutID = 1
                    Nmax = round(N*exp(Rmax * SM.tend))
                else
                    return SkinStemCellResults2D(scpool, SM)
                end
            end
        end
    end
    return SkinStemCellResults2D(scpool, SM)
end
