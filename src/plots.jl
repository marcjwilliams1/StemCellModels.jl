@recipe function f(results::StemCellResults; cutoff = 100)
    n = results.clonesize[:n][1:cutoff]
    Pn = results.clonesize[:Cn][1:cutoff]
    Pntheory = results.clonesize[:Cntheory]
    Pntheoryselection = results.clonesize[:Cntheoryselection]

    @series begin
        seriestype --> :bar
        fillcolor --> :darkslategrey
        legend --> false
         linecolor --> :white
        n, Pn
    end

    @series begin
        seriestype --> :line
        fillcolor --> :darkred
        legend --> false
        linewidth --> 3
        n, Pntheory
    end

    @series begin
        seriestype --> :line
        fillcolor --> :darkgreen
        legend --> false
        linewidth --> 3
        n, Pntheoryselection
    end
end
