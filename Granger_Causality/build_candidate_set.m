function candidates = build_candidate_set(Xq, maxLags, excludeIdx)
    [T, M] = size(Xq);
    candidates = [];
    for m = 1:M
        if m == excludeIdx
            continue
        end
        for lag = 1:maxLags
            candidates = [candidates; m, lag];
        end
    end
end
