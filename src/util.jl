function gethitchikers(scs::StemCellResults)
    idx = findall(x -> x == 0, scs.hitchikers)
    idx2 = findall(x -> x == 1, scs.hitchikers)
    return scs.mutationsize_p[idx], scs.mutationsize_p[idx2]
end
