% === GŁÓWNY SKRYPT: Analiza Odporności na Ataki (Alteration Attacks) ===
clear; clc; close all;

fprintf('Rozpoczynanie analizy odporności na ataki...\n');
fprintf('WYMAGA: Image Processing Toolbox (dla psnr, ssim, imnoise).\n\n');

% --- Konfiguracja ---
algorithms = {'AES-256-CBC', 'Blowfish-CBC', 'Chaotyczny (MATLAB)'};
attack_types = {'crop', 'salt_pepper', 'gaussian', 'tamper'};
% Definiujemy poziomy intensywności dla każdego ataku
intensities_crop = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5];
intensities_salt = [0.005, 0.01, 0.02, 0.05, 0.1];
intensities_gauss = [0.001, 0.005, 0.01, 0.02, 0.05];
intensities_tamper = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5];
% Zbiór intensywności dla pętli
intensity_sets = {intensities_crop, intensities_salt, intensities_gauss, intensities_tamper};

num_algos = length(algorithms);
num_attacks = length(attack_types);

% Wczytaj obrazy
try
    img_orig = imread('cameraman.tif'); 
catch
    img_orig = uint8(randi([0 255], 256, 256));
end
try
    % Wczytaj obraz-intruza dla ataku 'tamper'
    img_tamper = imread('peppers.png'); 
    img_tamper = imresize(img_tamper, size(img_orig));
    if size(img_tamper, 3) > 1
        img_tamper = rgb2gray(img_tamper);
    end
catch
    warning('Nie można wczytać "peppers.png". Generowanie losowego intruza.');
    img_tamper = uint8(randi([0 255], size(img_orig)));
end

img_flat = img_orig(:)';
[rows, cols] = size(img_orig);

% --- Generowanie kluczy (stałe dla wszystkich testów) ---
fprintf('Generowanie kluczy...\n');
% Klucze dla algorytmów klasycznych
key_aes = uint8(randi([0 255], 1, 32));
iv_aes = uint8(randi([0 255], 1, 16));
key_bf = uint8(randi([0 255], 1, 16));
iv_bf = uint8(randi([0 255], 1, 8));
% Klucze dla algorytmu chaotycznego (są zależne od obrazu, ale funkcja szyfrująca zwróci je)

% --- Szyfrowanie (raz, na początku) ---
fprintf('Szyfrowanie obrazu oryginalnego każdym algorytmem...\n');
% POPRAWKA: Musimy przyciąć dopełnienie (padding)
ct_aes_flat = crypto_mex(img_flat, key_aes, iv_aes, 'aes-256-cbc', true);
C_aes = reshape(ct_aes_flat(1:numel(img_orig)), [rows, cols]);

% POPRAWKA: Musimy przyciąć dopełnienie (padding)
ct_bf_flat = crypto_mex(img_flat, key_bf, iv_bf, 'bf-cbc', true);
C_bf = reshape(ct_bf_flat(1:numel(img_orig)), [rows, cols]);

[C_chaos, keys_chaos] = encrypt_circular_chaotic(img_orig);
% Zapisz szyfrogramy i klucze w komórkach
ciphertexts = {C_aes, C_bf, C_chaos};
keys = { {key_aes, iv_aes}, {key_bf, iv_bf}, keys_chaos };

% --- Pre-alokacja wyników ---
% Używamy komórek, bo wektory intensywności mają różne długości
results_psnr = cell(num_algos, num_attacks);
results_ssim = cell(num_algos, num_attacks);
results_sd = cell(num_algos, num_attacks);

% =========================================================================
% === GŁÓWNA PĘTLA TESTOWA (Algo -> Atak -> Intensywność) ===
% =========================================================================
for a_idx = 1:num_algos
    fprintf('\nTestowanie algorytmu: %s\n', algorithms{a_idx});
    
    for t_idx = 1:num_attacks
        attack = attack_types{t_idx};
        intensities = intensity_sets{t_idx};
        num_intensities = length(intensities);
        
        fprintf('  -> Atak: %s (%d poziomów intensywności)\n', attack, num_intensities);
        
        % Wektory na wyniki dla danego ataku
        psnr_vec = zeros(1, num_intensities);
        ssim_vec = zeros(1, num_intensities);
        sd_vec = zeros(1, num_intensities);
        
        C_orig = ciphertexts{a_idx}; % Weź oryginalny, czysty szyfrogram
        
        for i_idx = 1:num_intensities
            intensity = intensities(i_idx);
            
            % 1. Zastosuj atak
            C_attacked = apply_attack(C_orig, attack, intensity, img_tamper);
            
            % 2. Odszyfruj uszkodzony szyfrogram (z oryginalnym kluczem)
            img_decrypted = [];
            switch a_idx
                case 1 % AES
                    k = keys{a_idx}{1}; v = keys{a_idx}{2};
                    D_flat = crypto_mex(C_attacked(:)', k, v, 'aes-256-cbc', false);
                    img_decrypted = reshape(D_flat(1:numel(img_orig)), [rows, cols]);
                case 2 % Blowfish
                    k = keys{a_idx}{1}; v = keys{a_idx}{2};
                    D_flat = crypto_mex(C_attacked(:)', k, v, 'bf-cbc', false);
                    img_decrypted = reshape(D_flat(1:numel(img_orig)), [rows, cols]);
                case 3 % Chaotic
                    k_struct = keys{a_idx};
                    img_decrypted = decrypt_circular_chaotic(C_attacked, k_struct);
            end
            
            % 3. Oblicz metryki (porównując z oryginałem A)
            [psnr_val, ssim_val, sd_val] = calculate_alteration_metrics(img_orig, img_decrypted);
            psnr_vec(i_idx) = psnr_val;
            ssim_vec(i_idx) = ssim_val;
            sd_vec(i_idx) = sd_val;
        end
        
        % Zapisz wektory wyników do głównej komórki
        results_psnr{a_idx, t_idx} = psnr_vec;
        results_ssim{a_idx, t_idx} = ssim_vec;
        results_sd{a_idx, t_idx} = sd_vec;
    end
end
fprintf('\nAnaliza zakończona. Generowanie wykresów...\n');

% =========================================================================
% === WIZUALIZACJA WYNIKÓW (jak Fig. 24, 25) ===
% =========================================================================
figure;
set(gcf, 'Name', 'Analiza Odporności na Ataki', 'NumberTitle', 'off', 'WindowState', 'maximized');
marker_styles = {'-o', '-s', '-^'}; % Style dla AES, Blowfish, Chaotic

for t_idx = 1:num_attacks
    attack = attack_types{t_idx};
    intensities = intensity_sets{t_idx};
    
    % --- Wykres PSNR (Im wyżej, tym lepiej) ---
    subplot(3, num_attacks, t_idx);
    hold on;
    for a_idx = 1:num_algos
        plot(intensities, results_psnr{a_idx, t_idx}, marker_styles{a_idx}, 'LineWidth', 1.5);
    end
    hold off;
    title(sprintf('Atak: %s', strrep(attack, '_', ' ')));
    ylabel('PSNR (dB)');
    if t_idx == 1, legend(algorithms, 'Location', 'southwest'); end
    grid on;

    % --- Wykres SSIM (Im wyżej, tym lepiej) ---
    subplot(3, num_attacks, t_idx + num_attacks);
    hold on;
    for a_idx = 1:num_algos
        plot(intensities, results_ssim{a_idx, t_idx}, marker_styles{a_idx}, 'LineWidth', 1.5);
    end
    hold off;
    ylabel('SSIM');
    grid on;

    % --- Wykres SD (Im niżej, tym lepiej) ---
    subplot(3, num_attacks, t_idx + 2*num_attacks);
    hold on;
    for a_idx = 1:num_algos
        plot(intensities, results_sd{a_idx, t_idx}, marker_styles{a_idx}, 'LineWidth', 1.5);
    end
    hold off;
    ylabel('Spectral Distortion (SD)');
    xlabel('Intensywność ataku');
    grid on;
end

disp('Gotowe.');