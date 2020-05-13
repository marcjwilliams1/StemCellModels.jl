# StemCellModels.jl

julia package to simulate clonal expansions in normal tissues, the simulation models stem cells undergoing (biased) drift. This was used in the following paper:
[Measuring the distribution of fitness effects in somatic evolution by combining clonal dynamics with dN/dS ratios](https://elifesciences.org/articles/48714).

## Install package

```julia
] add StemCellModels
```

## Getting started
The first step is to setup a `StemCellModel` object which sets all the parameters of the simulation.
```julia
Ncells = 10^5
SM = StemCellModel(Ncells,
                   Δmut = 0.0, #"bias" of new driver mutations (ie selection coefficient)
                   μp = 0.1, #passenger mutation rate per effective division
                   μd = 0.1, # driver mutation rate per effective division
                   tend = 80.0, # time in years to simulate for
                   r = 0.5, #proportion of stem cells that result in symmetric division
                   λ = 1 # rate of stem cell division (per year)
                   )
```
Note that the effective of rate of division is given by `rλ`, so it's important to consider this product. I would suggest fixing one of these numbers and varying the other. For example `r = 0.5, λ = 1` will result in the same dynamics as `r = 1, λ = 0.5`. This is because we only consider stem cells in this model and do not include differentiated cells.

Then we can simulate a population of stem cells using the `runsimulation` functions.

```julia
sim = runsimulation(SM,
                    progress = true, #include progress bar
                    onedriver = true, #if this is set to true (default) then cells can only accumulate a single driver mutations
                    restart = true # restart simulation if population dies out
                    )
```
Note on the `onedriver` option, setting this to `true` means cells can only accumulate a single driver mutation. This is the regime in which the theoretical model described in the paper is valid. This assumption also seems to be mostly born out in the data but may be violated in some cases. If this is set to `false` then cells will gain multiple driver mutations but their fitness will be the same as if they had a single driver mutation.

You can access the output using the fields in results object. calling `sim` will also print out a summary.

```julia
sim.mutationfrequencies_d # frequencies of driver mutations
sim.mutationfrequencies_p # frequencies of passenger mutations

sim.mutationsize_d #number of cells each mutations is present int
sim.mutationsize_p

sim.hitchikers #binary array denoting whether each mutation is a hitchiker or not

```
