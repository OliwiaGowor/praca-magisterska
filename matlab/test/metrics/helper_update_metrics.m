function res = helper_update_metrics(res, alg_idx, run_idx, t_enc, t_dec, C1, C2, PT)
    % Calculates metrics and updates the main results structure
    
    % Record times
    res.times_enc(alg_idx, run_idx) = t_enc;
    res.times_dec(alg_idx, run_idx) = t_dec;
    
    % NPCR / UACI (only if dimensions match)
    if numel(C1) == numel(C2)
        [n, u] = calculate_npcr_uaci(C1, C2);
        res.npcr(alg_idx, run_idx) = n;
        res.uaci(alg_idx, run_idx) = u;
    end
    
    % Image statistics
    res.entropy(alg_idx, run_idx) = calculate_entropy(C1);
    [ch, cv, cd] = calculate_correlation(C1);
    res.corr_h(alg_idx, run_idx) = ch;
    res.corr_v(alg_idx, run_idx) = cv;
    res.corr_d(alg_idx, run_idx) = cd;
    
    % Save images (overwrites previous, keeps the last run)
    res.images_enc{alg_idx} = C1;
    res.images_dec{alg_idx} = PT;
end
