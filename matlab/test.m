% === GŁÓWNY SKRYPT PORÓWNAWCZY (Wersja 3.0 - z ChaCha20) ===
clear; clc; close all;

% --- Konfiguracja testów ---
N_RUNS = 10;
fprintf('Rozpoczynanie %d przebiegów testowych...\n', N_RUNS);

% --- Przygotowanie danych (obrazu) ---
try
    img_orig = imread('cameraman.tif'); 
catch
    img_orig = uint8(randi([0 255], 256, 256));
end
img_size = size(img_orig);

% --- Analiza obrazu oryginalnego ---
entropy_orig = calculate_entropy(img_orig);
[corr_h_orig, corr_v_orig, corr_d_orig] = calculate_correlation(img_orig);

% --- Obraz zmodyfikowany (dla NPCR/UACI) ---
img_mod = img_orig;
img_mod(1) = bitxor(img_mod(1), 1); % Zmień pierwszy piksel
img_flat = img_orig(:)';
img_flat_mod = img_mod(:)';

% --- Pre-alokacja macierzy na wyniki ---
% === NOWA LISTA ALGORYTMÓW ===
titles = {'AES-CBC (Hardware-C)', 'Blowfish-CBC (Software-C)', 'ChaCha20 (Software-C)', 'Chaotyczny (MATLAB/JIT)'};
num_algos = length(titles);

times_enc = NaN(num_algos, N_RUNS);
times_dec = NaN(num_algos, N_RUNS);
npcr_values = NaN(num_algos, N_RUNS);
uaci_values = NaN(num_algos, N_RUNS);
entropy_values = NaN(num_algos, N_RUNS);
corr_h_values = NaN(num_algos, N_RUNS);
corr_v_values = NaN(num_algos, N_RUNS);
corr_d_values = NaN(num_algos, N_RUNS);

images_enc = cell(num_algos, 1);
images_dec = cell(num_algos, 1);

% =========================================================================
% === GŁÓWNA PĘTLA TESTOWA ===
% =========================================================================
for run_idx = 1:N_RUNS
    fprintf('--- Rozpoczynanie przebiegu %d / %d ---\n', run_idx, N_RUNS);
    
    % === TEST 1: AES-256-CBC (C/MEX) ===
    key_aes = uint8(randi([0 255], 1, 32));
    iv_aes = uint8(randi([0 255], 1, 16)); 
    tic;
    ct_aes_C1_flat = crypto_mex(img_flat, key_aes, iv_aes, 'aes-256-cbc', true);
    times_enc(1, run_idx) = toc; 
    ct_aes_C2_flat = crypto_mex(img_flat_mod, key_aes, iv_aes, 'aes-256-cbc', true);
    tic;
    pt_aes_padded = crypto_mex(ct_aes_C1_flat, key_aes, iv_aes, 'aes-256-cbc', false);
    times_dec(1, run_idx) = toc;
    C1 = reshape(ct_aes_C1_flat(1:numel(img_orig)), img_size);
    C2 = reshape(ct_aes_C2_flat(1:numel(img_orig)), img_size);
    [npcr_values(1, run_idx), uaci_values(1, run_idx)] = calculate_npcr_uaci(C1, C2);
    entropy_values(1, run_idx) = calculate_entropy(C1);
    [corr_h_values(1, run_idx), corr_v_values(1, run_idx), corr_d_values(1, run_idx)] = calculate_correlation(C1);
    images_enc{1} = C1; 
    images_dec{1} = reshape(pt_aes_padded(1:numel(img_orig)), img_size);

    % === TEST 2: Blowfish-CBC (C/MEX) ===
    key_bf = uint8(randi([0 255], 1, 16));
    iv_bf = uint8(randi([0 255], 1, 8)); 
    tic;
    ct_bf_C1_flat = crypto_mex(img_flat, key_bf, iv_bf, 'bf-cbc', true);
    times_enc(2, run_idx) = toc;
    ct_bf_C2_flat = crypto_mex(img_flat_mod, key_bf, iv_bf, 'bf-cbc', true);
    tic;
    pt_bf_padded = crypto_mex(ct_bf_C1_flat, key_bf, iv_bf, 'bf-cbc', false);
    times_dec(2, run_idx) = toc;
    C1 = reshape(ct_bf_C1_flat(1:numel(img_orig)), img_size);
    C2 = reshape(ct_bf_C2_flat(1:numel(img_orig)), img_size);
    [npcr_values(2, run_idx), uaci_values(2, run_idx)] = calculate_npcr_uaci(C1, C2);
    entropy_values(2, run_idx) = calculate_entropy(C1);
    [corr_h_values(2, run_idx), corr_v_values(2, run_idx), corr_d_values(2, run_idx)] = calculate_correlation(C1);
    images_enc{2} = C1; 
    images_dec{2} = reshape(pt_bf_padded(1:numel(img_orig)), img_size);

    % === TEST 3: ChaCha20 (C/MEX) ===
    key_chacha = uint8(randi([0 255], 1, 32)); % 32 bajty
    iv_chacha = uint8(randi([0 255], 1, 16));  % OpenSSL EVP używa 16-bajtowego IV
    tic;
    ct_chacha_C1_flat = crypto_mex(img_flat, key_chacha, iv_chacha, 'chacha20', true);
    times_enc(3, run_idx) = toc;
    ct_chacha_C2_flat = crypto_mex(img_flat_mod, key_chacha, iv_chacha, 'chacha20', true);
    tic;
    pt_chacha_flat = crypto_mex(ct_chacha_C1_flat, key_chacha, iv_chacha, 'chacha20', true); % Szyfrowanie to deszyfrowanie
    times_dec(3, run_idx) = toc;
    C1 = reshape(ct_chacha_C1_flat, img_size);
    C2 = reshape(ct_chacha_C2_flat, img_size);
    [npcr_values(3, run_idx), uaci_values(3, run_idx)] = calculate_npcr_uaci(C1, C2);
    entropy_values(3, run_idx) = calculate_entropy(C1);
    [corr_h_values(3, run_idx), corr_v_values(3, run_idx), corr_d_values(3, run_idx)] = calculate_correlation(C1);
    images_enc{3} = C1; 
    images_dec{3} = reshape(pt_chacha_flat, img_size);
    
    % === TEST 4: ALGORYTM CHAOTYCZNY (MATLAB .m) ===
    try
        tic;
        [C1_chaos, keys_chaos_C1] = encrypt_circular_chaotic(img_orig);
        times_enc(4, run_idx) = toc;
        [C2_chaos, keys_chaos_C2] = encrypt_circular_chaotic(img_mod);
        tic;
        pt_chaos = decrypt_circular_chaotic(C1_chaos, keys_chaos_C1);
        times_dec(4, run_idx) = toc;
        [npcr_values(4, run_idx), uaci_values(4, run_idx)] = calculate_npcr_uaci(C1_chaos, C2_chaos);
        entropy_values(4, run_idx) = calculate_entropy(C1_chaos);
        [corr_h_values(4, run_idx), corr_v_values(4, run_idx), corr_d_values(4, run_idx)] = calculate_correlation(C1_chaos);
        images_enc{4} = C1_chaos; 
        images_dec{4} = pt_chaos;
    catch ME
        warning('Przebieg %d dla algorytmu chaotycznego nie powiódł się: %s', run_idx, ME.message);
    end
end
fprintf('--- Wszystkie przebiegi zakończone! ---\n\n');

% =========================================================================
% === OBLICZANIE ŚREDNICH I ZAPIS DO PLIKU ===
% =========================================================================

avg_times_enc = mean(times_enc, 2, 'omitnan');
avg_times_dec = mean(times_dec, 2, 'omitnan');
avg_npcr = mean(npcr_values, 2, 'omitnan');
avg_uaci = mean(uaci_values, 2, 'omitnan');
avg_entropy = mean(entropy_values, 2, 'omitnan');
avg_corr_h = mean(corr_h_values, 2, 'omitnan');
avg_corr_v = mean(corr_v_values, 2, 'omitnan');
avg_corr_d = mean(corr_d_values, 2, 'omitnan');

filename = 'encryption_results_MATLAB.mat';
save(filename, 'titles', 'N_RUNS', ...
     'times_enc', 'times_dec', ...
     'npcr_values', 'uaci_values', ...
     'entropy_values', ...
     'corr_h_values', 'corr_v_values', 'corr_d_values', ...
     'entropy_orig', 'corr_h_orig', 'corr_v_orig', 'corr_d_orig', ...
     'avg_times_enc', 'avg_times_dec', 'avg_npcr', 'avg_uaci', ...
     'avg_entropy', 'avg_corr_h', 'avg_corr_v', 'avg_corr_d');
fprintf('Wszystkie surowe wyniki (%d przebiegów) zostały zapisane w pliku: %s\n', N_RUNS, filename);
% =========================================================================
% === WIZUALIZACJA ŚREDNICH WYNIKÓW (POPRAWIONA SKŁADNIA) ===
% =========================================================================
disp('Generowanie wizualizacji na podstawie średnich wyników...');
figure;
set(gcf, 'Name', 'Porównanie: Wydajność (Wykres) i Bezpieczeństwo (Obrazy)', 'NumberTitle', 'off', 'WindowState', 'maximized');

subplot_height = num_algos + 2; % Teraz 4 + 2 = 6

% --- Górny panel: Wykres wydajności (ŚREDNIE) ---
subplot(subplot_height, 3, [2 3]);
bar_data = [avg_times_enc avg_times_dec];
b = bar(categorical(titles), bar_data);
legend('Szyfrowanie', 'Deszyfrowanie', 'Location', 'northwest');
title(sprintf('Porównanie WYDAJNOŚCI (Średnia z %d przebiegów, Obraz: %dx%d)', N_RUNS, img_size(1), img_size(2)));
ylabel('Czas (sekundy)');
if ~isempty(avg_times_enc) && max(avg_times_enc) / min(avg_times_enc(avg_times_enc>0)) > 100
    set(gca, 'YScale', 'log'); 
    ylabel('Czas (sekundy) - Skala Log');
end
grid on;

% --- Dolne panele: Obrazy ---
subplot(subplot_height, 3, 4);
imshow(img_orig);
title_str = sprintf('Oryginał\nEntropia: %.4f\nKorelacja H/V/D: %.3f / %.3f / %.3f', ...
                    entropy_orig, corr_h_orig, corr_v_orig, corr_d_orig);
title(title_str);
ylabel('WEJŚCIE');

subplot(subplot_height, 3, 7);
imshow(img_mod);
title('Oryginał (zmieniony 1px)');
ylabel('WEJŚCIE MOD');

% Wykresy dla każdego algorytmu
for i = 1:num_algos
    subplot(subplot_height, 3, (i+1)*3 + 2); % Poprawione indeksowanie
    imshow(images_enc{i});
    title_str = sprintf('Szyfrogram (Średni Enc: %.4f s)\nŚr. Ent: %.4f, Śr. NPCR: %.2f%%, Śr. UACI: %.2f%%\nŚr. Korelacja H/V/D: %.3f / %.3f / %.3f', ...
                         avg_times_enc(i), ...
                         avg_entropy(i), avg_npcr(i), avg_uaci(i), ...
                         avg_corr_h(i), avg_corr_v(i), avg_corr_d(i));
    title(title_str, 'FontSize', 8); % Mniejsza czcionka
    ylabel(titles{i}, 'FontWeight', 'bold');
    
    subplot(subplot_height, 3, (i+1)*3 + 3); % Poprawione indeksowanie
    imshow(images_dec{i});
    title(sprintf('Deszyfrogram\n(Średni Dec: %.4f s)', avg_times_dec(i)));
end

% === POPRAWKA BŁĘDU ===
% Ukryj puste sloty używając oddzielnych poleceń
subplot(subplot_height, 3, 1); axis('off'); % Zamiast .axis('off')
subplot(subplot_height, 3, 2); axis('off'); % (puste miejsce nad wykresem)
subplot(subplot_height, 3, 3); axis('off'); % (puste miejsce nad wykresem)
subplot(subplot_height, 3, 10); axis('off'); % (puste miejsce obok AES)
subplot(subplot_height, 3, 13); axis('off'); % (puste miejsce obok Blowfish)
% Dodano brakujące ukrywanie dla pozostałych algorytmów
subplot(subplot_height, 3, 16); axis('off'); % (puste miejsce obok ChaCha20)

disp('Gotowe.');