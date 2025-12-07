function avg = process_results(raw, data, config)
    % Calculates averages and saves results to a .mat file
    
    avg.times_enc = mean(raw.times_enc, 2, 'omitnan');
    avg.times_dec = mean(raw.times_dec, 2, 'omitnan');
    avg.npcr = mean(raw.npcr, 2, 'omitnan');
    avg.uaci = mean(raw.uaci, 2, 'omitnan');
    avg.entropy = mean(raw.entropy, 2, 'omitnan');
    avg.corr_h = mean(raw.corr_h, 2, 'omitnan');
    avg.corr_v = mean(raw.corr_v, 2, 'omitnan');
    avg.corr_d = mean(raw.corr_d, 2, 'omitnan');

    filename = sprintf('results/encryption_results_MATLAB_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
    
    % Prepare local variables for saving
    titles = config.titles;
    N_RUNS = config.N_RUNS;
    
    % Save variables to disk
    save(filename, 'titles', 'N_RUNS', 'raw', 'avg', 'data');
    fprintf('Wyniki zapisano w: %s\n', filename);
end
