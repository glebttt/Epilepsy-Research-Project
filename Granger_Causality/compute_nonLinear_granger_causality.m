%% NON-Linear granger causality function

function results = compute_nonLinear_granger_causality(X, L, Q, nSurr, alpha)

    [T, M] = size(X);
    % Normalize to zero mean, unit variance
    Xn = zscore(X);
    Xq = quantize_uniform(Xn, Q);   
    Cmat   = zeros(M,M);
    pMat   = ones(M,M);
    sigMat = zeros(M,M);
    Kfull  = zeros(M,M);
    Kred   = zeros(M,M);

    for i = 1:M
        for j = 1:M
            
            if i == j
                continue;
            end
            
            %% ---- Full candidate set Ω ----
            cand_full = build_candidate_set(Xq, L, []);      % all variables
            cand_red  = build_candidate_set(Xq, L, j);       % excluding Xj
            
            %% ---- Nonuniform embedding ----
            [Kvec_full, Hfull] = nonuniform_embedding_entropy(Xq(:,i), cand_full, Q);
            [Kvec_red, Hred]   = nonuniform_embedding_entropy(Xq(:,i), cand_red, Q);
            
            Kfull(i,j) = length(Kvec_full);
            Kred(i,j)  = length(Kvec_red);
            
            %% ---- Corrected causality ----
            C_real = 1 - Hfull / Hred;
            Cmat(i,j) = C_real;
            
            %% ========== Surrogate testing ==========
            C_surr = zeros(nSurr,1);
            
            for s = 1:nSurr
                
                Xq_surr = Xq;
                Xq_surr(:,j) = time_shift_surrogate(Xq(:,j));
                
                cand_full_s = build_candidate_set(Xq_surr, L, []);
                
                [~, Hfull_s] = nonuniform_embedding_entropy(Xq_surr(:,i), cand_full_s, Q);
                                
                C_surr(s) = 1 - Hfull_s / Hred; % reduced stays original
            end
            
            %% ---- Rank test ----
            allC = [C_real; C_surr];
            [~, idx] = sort(allC);
            rank = find(idx == 1);   % position of real value
            
            pval = 1 - (rank - 0.326)/(nSurr + 1 + 0.348);
            
            pMat(i,j)   = pval;
            sigMat(i,j) = pval < alpha;
            
        end
    end
    results.C      = Cmat;
    results.p      = pMat;
    results.sig    = sigMat;
    results.Kfull  = Kfull;
    results.Kred   = Kred;
end