using Test
using StemCellModels
using Statistics
using Random

tests = ["skin"]

println("Running tests ...")

@time for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
