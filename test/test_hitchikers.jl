Random.seed!(123)

Nsims = 10^2
SM = SkinStemCellModel(100, Δmut = 0.0, μp = 0.05, μd = 0.05, tend = 30.0, r = 1.0, λ = 0.5)

global neutral = Float64[]
global hitchikers = Float64[]
global drivers = Float64[]
for i in 1:Nsims
    scst = runsimulation(SM)
    x = gethitchikers(scst)
    append!(drivers, scst.mutationsize_d)
    append!(neutral, x[1])
    append!(hitchikers, x[2])
end


SM = SkinStemCellModel(100, Δmut = 0.2, μp = 0.05, μd = 0.05, tend = 30.0, r = 1.0, λ = 0.5)

global neutral2 = Float64[]
global hitchikers2 = Float64[]
global drivers2 = Float64[]
for i in 1:Nsims
    scst = runsimulation(SM)
    x = gethitchikers(scst)
    append!(drivers2, scst.mutationsize_d)
    append!(neutral2, x[1])
    append!(hitchikers2, x[2])
end

Random.seed!(123)
SM = SkinStemCellModel(100, Δmut = 0.2, μp = 0.2, μd = 0.2, tend = 30.0, r = 1.0, λ = 0.5)
scst = runsimulation(SM)
x = gethitchikers(scst)


@test mean(hitchikers2) > mean(neutral2)
@test mean(drivers2) > mean(hitchikers2)
@test mean(drivers2) > mean(neutral2)
@test mean(drivers2) > mean(drivers)
@test mean(hitchikers2) > mean(hitchikers)
