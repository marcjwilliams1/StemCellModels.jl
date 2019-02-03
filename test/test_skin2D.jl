Random.seed!(1)

Ns = 100
SM = SkinStemCellModel2D(Ns)

pool = StemCellModels.createstemcellpool((SM))
# x,y coordinates of cells should all be unique
@test length(unique(keys(pool))) == Ns

(newx, newy) = StemCellModels.getfreespace(pool, (50, 50))

x = runsimulation(SM)
