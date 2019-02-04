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

    function SkinStemCellResults2D(stemcells, SM)
        idx = values(stemcells) .!= 0
        return new(collect(values(stemcells))[idx], SM)
    end
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
    stemcell = SkinStemCell2D(SM.r, SM.λ, SM.Δ, SM.μd, SM.μp, Int64[], Int64[], midpoint[1], midpoint[2])
    stemcellpool = DefaultDict{Tuple, Any}(0)
    push!(stemcellpool, (midpoint[1], midpoint[2]) => stemcell)
    Nstemcells = 1
    while Nstemcells < SM.Ns
        r = radius * sqrt(rand())
        θ = rand() * 2 * π
        x = round(midpoint[1] + r * cos(θ))
        y = round(midpoint[2] + r * sin(θ))
        if stemcellpool[(x,y)] == 0
            push!(stemcellpool,
            (x, y) => SkinStemCell2D(SM.r, SM.λ, SM.Δ, SM.μd, SM.μp, Int64[], Int64[], x, y))
            Nstemcells += 1
        end
    end
    return stemcellpool
end

function getfreespace(scpool, randomsc)
    space = Tuple[]
    for x in [-1, 1]
        for y in [1, 1]
            newx = scpool[randomsc].x + x
            newy = scpool[randomsc].y + y
            if scpool[(newx,newy)] == 0
                push!(space, (newx, newy))
            end
        end
    end

    if isempty(space)
        direction = rand([-1, 1], 2)
        newx = scpool[randomsc].x + direction[1]
        newy = scpool[randomsc].y + direction[2]
        return (newx, newy)
    else
        return rand(space)
    end
end

function placestemcell(scpool::DefaultDict, randomsc)
    (newx, newy) = getfreespace(scpool, randomsc)
    #println("Current stem cell: ($(scpool[randomsc].x), $(scpool[randomsc].y), index:$(randomsc))")
    #println("New stem cell: ($(newx), $(newy))")
    #println("Size of scpool: $(length(scpool))")
    if scpool[(newx,newy)] == 0
        newstemcell = copyskinstemcell(scpool[randomsc])
        newstemcell.x = newx
        newstemcell.y = newy
        push!(scpool, (newx, newy) => newstemcell)
        return scpool, (newx, newy)
    else
        delete!(scpool, (newx, newy))
        newstemcell = copyskinstemcell(scpool[randomsc])
        newstemcell.x = newx
        newstemcell.y = newy
        push!(scpool, (newx, newy) => newstemcell)
        return scpool, (newx, newy)
    end
end

function cellturnover(scpool::DefaultDict, randomsc, mutID, SM; onedriver = true)
    r = rand()
    #println("Number of cells preturnover: $(length(scpool))")
    if r < ((1 + scpool[randomsc].Δ) / 2)
        scpool, index = placestemcell(scpool, randomsc)
        scpool[index], mutID =
            newmutations(scpool[index], mutID, SM, onedriver = onedriver)
    else
        delete!(scpool, randomsc)
    end
    #println("Number of cells posturnover: $(length(scpool))")
    return scpool, mutID
end

function randomstemcell(scpool)
    x = 0
    while x == 0
        rsc = rand(keys(scpool))
        #println(scpool)
        println(length(keys(scpool)))
        if scpool[rsc] != 0
            x = 1
            return rsc
        end
    end
end

function runsimulation(SM::SkinStemCellModel2D; progress = false, restart = false, finish = "time", onedriver = true)
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
            randomsc = randomstemcell(scpool)
            u1 = rand(Uniform(0, Rmax))
            if u1 < (scpool[randomsc].λ * scpool[randomsc].r * 2)
                #println("divide")
                scpool, mutID = cellturnover(scpool, randomsc,
                                        mutID, SM,
                                        onedriver = onedriver)
                Δt = (1 / (Rmax * N)) * log(1/rand())
                t += Δt
                N = sum(values(scpool) .!= 0)
                Ndivs += 1
            else
                #println("die")
                Δt = (1 / (Rmax * N)) * log(1/rand())
                t += Δt
                N = sum(values(scpool) .!= 0)
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
                    return SkinStemCellResults2D(scpool, SM)
                end
            end
        end
    else
        #println("Nmax is $(Nmax)")
        while N < Nmax
            randomsc = rand(keys(scpool))
            u1 = rand(Uniform(0, Rmax))
            if u1 < (scpool[randomsc].λ * scpool[randomsc].r * 2)
                #println("divide")
                scpool, mutID = cellturnover(scpool, randomsc,
                                        mutID, SM,
                                        onedriver = onedriver)
                Δt = (1 / (Rmax * N)) * log(1/rand())
                t += Δt
                N = sum(values(pool) .!= 0)
                Ndivs += 1
            else
                #println("die")
                Δt = (1 / (Rmax * N)) * log(1/rand())
                t += Δt
                N = sum(values(pool) .!= 0)
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
                    return SkinStemCellResults2D(scpool, SM)
                end
            end
        end
    end
    return SkinStemCellResults2D(scpool, SM)
end
