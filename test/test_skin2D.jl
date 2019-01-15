Random.seed!(1)

Ns = 100
SM = SkinStemCellModel2D(Ns)

pool, scids = StemCellModels.createstemcellpool((SM))
# x,y coordinates of cells should all be unique
@test length(unique(keys(scids))) == Ns

runsimulation(SM)
