function out_flat = hc_confusion_inv(in_flat, seq_conf)
    % HC_CONFUSION_INV Odwrócenie permutacji
    
    [~, sort_idx] = sort(seq_conf);
    
    ranks = zeros(size(seq_conf));
    ranks(sort_idx) = 1:length(seq_conf);
    
    out_flat = zeros(size(in_flat));
    
    for i = 1:length(in_flat)
        idx = ranks(i);
        out_flat(idx) = in_flat(i);
    end
end