using StemCellModels
SM = SkinStemCellModel2D(1000)
pool, scids = StemCellModels.createstemcellpool((SM))
x = map(x -> x.x, pool)
y = map(x -> x.y, pool)
[x y]
using DataFrames
using RCall
DF = DataFrame(x = x, y = y, z = 1)
@rput DF
R"library(ggplot2)"
R"ggplot(DF, aes(x = x, y = y, fill = z)) + geom_tile()"
