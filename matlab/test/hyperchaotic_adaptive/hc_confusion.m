function out_flat = hc_confusion(in_flat, seq_conf)
    % HC_CONFUSION Permutacja oparta na sortowaniu sekwencji chaotycznej
    
    [~, sort_idx] = sort(seq_conf);
    
    % Rangi elementów
    ranks = zeros(size(seq_conf));
    ranks(sort_idx) = 1:length(seq_conf);
    
    out_flat = zeros(size(in_flat));
    
    for i = 1:length(in_flat)
        idx = ranks(i); 
        out_flat(i) = in_flat(idx);
    end
end