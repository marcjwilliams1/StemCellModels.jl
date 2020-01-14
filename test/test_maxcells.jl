Random.seed!(123)

Nsims = 10^2
SM = StemCellModel(100, Δmut = 0.2, μp = 0.05, μd = 0.05, tend = 80.0, r = 1.0, λ = 0.5)
@time scst = runsimulation(SM, maxcells = 10^6);
