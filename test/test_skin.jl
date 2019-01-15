#test create stem cell pool
Random.seed!(1)

SM = SkinStemCellModel(100)

@test length(StemCellModels.createstemcellpool(SM)) == 100

# when Δ = 0.0 average population should be equal to starting population
N = 10
SM = SkinStemCellModel(N, Δ = 0.0)
popsize = Int64[]
for i in 1:10^4
    scs = runsimulation(SM)
    push!(popsize, length(scs.stemcells))
end
@test isapprox(round(mean(popsize)), Float64(N))

#variance should be equal to 2 * λ * r * N
@test isapprox(var(popsize), SM.r * SM.λ * SM.tend * 2 * N, rtol = 0.1)

#test probability distribution matches theory
N = 10^4
SM = SkinStemCellModel(N, Δ = 0.0, μp = 0.1, μd = 0.0, tend = 200.0)
results = runsimulation(SM, progress = true)
@test isapprox(round.(results.clonesize[:Pn], digits = 4),
 round.(results.clonesize[:Pntheory], digits = 4), rtol = 0.1)

#look at positive selection,
#average populaiton size should be N*exp(2rλΔt)
# and variance = (N / Δ) * exp(2rλΔt) * (exp(2rλΔt) - 1)
N = 10
SM = SkinStemCellModel(N, Δ = 0.1, tend = 10.0)
popsize = Int64[]
for i in 1:10^4
    scs = runsimulation(SM)
    push!(popsize, length(scs.stemcells))
end
Nexp = N * exp(2 * SM.r * SM.λ * SM.Δ * SM.tend)
@test isapprox(round(mean(popsize)), Nexp, rtol = 0.1)

varexp = (N / SM.Δ) * exp(2 * SM.r * SM.λ * SM.Δ * SM.tend) *
    (exp(2 * SM.r * SM.λ * SM.Δ * SM.tend) - 1)
@test isapprox(var(popsize), varexp, rtol = 0.1)

#test probability distribution for non-netural mutations

N = 10^4
SM = SkinStemCellModel(N, Δ = 0.01, μp = 0.1, μd = 0.0, tend = 200.0)
results = runsimulation(SM, progress = true)
isapprox(round.(results.clonesize[:Pn], digits = 4),
 round.(results.clonesize[:Pntheory], digits = 4), rtol = 0.1)

N = 10^4
SM = SkinStemCellModel(N, Δ = 0.0, Δmut = 0.1, μp = 0.0, μd = 0.1, tend = 50.0)
results = runsimulation(SM, progress = true, onedriver = true)
results.clonesize
isapprox(round.(results.clonesize[:Pn], digits = 4),round.(results.clonesize[:Pntheory], digits = 4), rtol = 0.1)

N = 10^3
SM = SkinStemCellModel(N, Δ = 0.0, μp = 1.0, μd = 0.0, tend = 200.0)
resultsneut = runsimulation(SM, progress = true)

SM = SkinStemCellModel(N, Δ = 0.1, μp = 1.0, μd = 0.0, tend = 2.0)
println(2 * SM.r * SM.λ * SM.Δ * SM.tend)
resultssel = runsimulation(SM, progress = true)


popsize = Int64[]
N = 1
SM = SkinStemCellModel(N, Δ = 0.1, μp = 0.0, μd = 0.0, tend = 20.0)
for i in 1:10^3
    resultssel = runsimulation(SM, progress = false)
    push!(popsize, length(resultssel.stemcells))
end
Nexp = N * exp(2 * SM.r * SM.λ * SM.Δ * SM.tend)


# Test when s != 0
# when Δ = 0.0 average population should be equal to starting population
N = 100
SMneutral = SkinStemCellModel(N, Δ = 0.0, Δmut = 0.0, μp = 0.0, μd = 0.0, s = 0.5)
SMselection = SkinStemCellModel(N, Δ = 0.0, Δmut = 0.01, μp = 0.0, μd = 1.0, s = 0.5)
SMselection2 = SkinStemCellModel(N, Δ = 0.0, Δmut = 0.1, μp = 0.0, μd = 1.0, s = 0.5)
SMselection3= SkinStemCellModel(N, Δ = 0.0, Δmut = 0.1, μp = 0.0, μd = 1.0, s = 1.0)
popsizeneutral = Int64[]
popsizeselection = Int64[]
popsizeselection2 = Int64[]
popsizeselection3 = Int64[]
divsneutral = Int64[]
divsselection = Int64[]
divsselection2= Int64[]
divsselection3 = Int64[]
for i in 1:10^4
    scsneutral = runsimulation(SMneutral, restart = true)
    scsselection = runsimulation(SMselection, restart = true)
    scsselection2 = runsimulation(SMselection2, restart = true)
    scsselection3 = runsimulation(SMselection3, restart = true)
    push!(popsizeneutral, length(scsneutral.stemcells))
    push!(popsizeselection, length(scsselection.stemcells))
    push!(popsizeselection2, length(scsselection2.stemcells))
    push!(popsizeselection3, length(scsselection3.stemcells))
end
mean(popsizeneutral)
mean(popsizeselection)
mean(popsizeselection2)
mean(popsizeselection3)
@test mean(popsizeselection3) > mean(popsizeselection2) > mean(popsizeselection) > mean(popsizeneutral)

SMneutral.λ * SMneutral.r * SMselection.tend  + 1
SMselection.λ * SMselection.r * (SMselection.s + 1) * SMselection.tend  + 1

N = 100
SMneutral = SkinStemCellModel(N, Δ = 0.0, μp = 0.0, μd = 0.0, s = 0.5, tend = 100.0)
SMselection = SkinStemCellModel(N, Δ = 0.0, μp = 0.0, μd = 0.1, s = 0.5, tend = 100.0)
scsneutral = runsimulation(SMneutral, restart = true)
scsselection = runsimulation(SMselection, restart = true)

length(scsneutral.stemcells)
length(scsselection.stemcells)
sum(map(x -> x.λ, scsselection.stemcells) .== 1.5)


# Test when s != 0
# when Δ = 0.0 average population should be equal to starting population
N = 3
SMneutral = SkinStemCellModel(N, Δ = 0.0, Δmut = 0.1, μp = 0.0, μd = 0.1, s = 0.0)
SMselection = SkinStemCellModel(N, Δ = 0.0, Δmut = 0.1, μp = 0.0, μd = 0.1, s = 1.0)
SMselection = SkinStemCellModel(N, Δ = 0.0, Δmut = 0.1, μp = 0.0, μd = 0.1, s = 1.0)
popsizeneutral = Int64[]
popsizeselection = Int64[]
divsneutral = Int64[]
divsselection = Int64[]
for i in 1:10^4
    scsneutral, ndneutral = runsimulation(SMneutral, restart = true)
    scsselection, ndselection = runsimulation(SMselection, restart = true)
    push!(popsizeneutral, length(scsneutral.stemcells))
    push!(popsizeselection, length(scsselection.stemcells))
    push!(divsneutral, ndneutral)
    push!(divsselection, ndselection)
end
mean(popsizeneutral)
mean(popsizeselection)
