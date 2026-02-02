% Nazwa skryptu: Jednowymiarowe_odwzorowanie_kosinusowo_wielomianowe.m
% Opis: Generuje diagram bifurkacyjny jako heatmapę oraz wykres wykładnika Lapunowa.

clear; clc; close all;

%% 1. Ustawienia Parametrów
% Parametry zgodne z sekcją 2 artykułu (Fig. 2 i Fig. 5)
mu_min = 0;
mu_max = 5;
n_points = 1000;            % Rozdzielczość parametru mu (oś X heatmapy)
mu_values = linspace(mu_min, mu_max, n_points);

a = 4;                      % Stała wartość parametru a
x0 = 0.1;                   % Warunek początkowy

iterations = 2000;          % Całkowita liczba iteracji na punkt
transient = 1000;           % Odrzucane iteracje początkowe (stan nieustalony)

% Ustawienia Heatmapy (rozdzielczość osi Y - wartości x)
y_res = 800;                % Liczba "koszyków" (bins) w pionie
y_limits = [-1, 1];         % Zakres wartości funkcji cosinus
y_edges = linspace(y_limits(1), y_limits(2), y_res + 1); % Krawędzie koszyków

% Inicjalizacja macierzy gęstości (Heatmapy) i wektora LE
heatmap_matrix = zeros(y_res, n_points);
lyapunov_exponents = zeros(1, n_points);

%% 2. Główna pętla obliczeniowa
fprintf('Generowanie danych dla heatmapy i LE (mu = [%.1f, %.1f])...\n', mu_min, mu_max);

for k = 1:n_points
    mu = mu_values(k);
    x = x0;
    
    % --- Iteracje przejściowe (dojście do atraktora) ---
    for i = 1:transient
        x = cos(mu * (x^3 + x) + a);
    end
    
    % --- Właściwe iteracje (zbieranie danych) ---
    current_orbit = zeros(1, iterations - transient);
    log_sum = 0; % Zmienna do sumowania logarytmów pochodnych (dla LE)
    
    for i = 1:(iterations - transient)
        x_prev = x;
        
        % Równanie mapy: x_k = cos(mu * (x^3 + x) + a)
        % Argument cosinusa oznaczmy jako u = mu*(x^3 + x) + a
        u_arg = mu * (x_prev^3 + x_prev) + a;
        x = cos(u_arg);
        
        current_orbit(i) = x;
        
        % Obliczenie pochodnej do Wykładnika Lapunowa
        % f'(x) = -sin(u) * u' = -sin(u) * mu * (3x^2 + 1)
        derivative = -sin(u_arg) * mu * (3 * x_prev^2 + 1);
        log_sum = log_sum + log(abs(derivative));
    end
    
    % --- 3. Budowanie Heatmapy ---
    % Obliczamy histogram wartości 'x' dla aktualnego 'mu'
    % Funkcja histcounts (dostępna w nowszych wersjach MATLABa) lub histc
    if exist('histcounts', 'file')
        counts = histcounts(current_orbit, y_edges);
    else
        % Kompatybilność ze starszymi wersjami
        counts = histc(current_orbit, y_edges);
        counts = counts(1:end-1); 
    end
    
    % Zapisujemy kolumnę do macierzy (transpozycja, aby pasowała do kolumny)
    heatmap_matrix(:, k) = counts';
    
    % --- 4. Obliczenie LE ---
    lyapunov_exponents(k) = log_sum / (iterations - transient);
end

%% 5. Rysowanie Wykresów

% --- Wykres 1: Diagram Bifurkacyjny (Heatmapa) ---
fig1 = figure('Name', 'Diagram Bifurkacyjny', 'Color', 'w', 'Position', [100, 100, 800, 500]);

% Używamy skali logarytmicznej (log1p), aby lepiej widzieć rzadkie punkty
imagesc(mu_values, y_limits, log1p(heatmap_matrix));

set(gca, 'YDir', 'normal'); % Odwrócenie osi Y (dół to -1, góra to 1)
colormap(jet);              % Mapa kolorów 'hot' (czarny tło, czerwone/żółte punkty)
colorbar;                   % Pasek legendy kolorów
caxis([0 max(max(log1p(heatmap_matrix)))]); % Skalowanie kolorów

xlabel(' Parametr kontrolny \mu');
ylabel('Zmienna stanu x');
title(['Diagram Bifurkacyjny - Jednowymiarowe odwzorowanie kosinusowo–wielomianowe (a = ' num2str(a) ')']);

% --- Wykres 2: Wykładnik Lapunowa ---
fig2 = figure('Name', 'Wykładnik Lapunowa', 'Color', 'w', 'Position', [150, 150, 800, 400]);

plot(mu_values, lyapunov_exponents, 'b-', 'LineWidth', 1.0);
hold on;
plot([mu_min, mu_max], [0, 0], 'r--', 'LineWidth', 1.5); % Linia graniczna chaosu (0)

xlim([mu_min mu_max]);
% Automatyczne dopasowanie osi Y
max_le = max(lyapunov_exponents);
min_le = min(lyapunov_exponents);
ylim([min_le - 0.1*(abs(min_le)), max_le + 0.1*(abs(max_le))]);

xlabel('\mu');
ylabel('Wykładnik Lapunowa (LE)');
title(['Wykładnik Lapunowa - Jednowymiarowe odwzorowanie kosinusowo–wielomianowe (a = ' num2str(a) ')']);
legend('LE', 'Granica chaosu (LE=0)', 'Location', 'best');
grid on;

%% 6. Zapisywanie do plików
filename_base = 'Jednowymiarowe_odwzorowanie_kosinusowo_wielomianowe';

% Zapis heatmapy
file_bif = [filename_base '_Bifurkacja_Heatmapa.png'];
saveas(fig1, file_bif);
fprintf('Zapisano heatmapę bifurkacyjną do: %s\n', file_bif);

% Zapis wykresu LE
file_lap = [filename_base '_Lapunov.png'];
saveas(fig2, file_lap);
fprintf('Zapisano wykres Lapunowa do: %s\n', file_lap);

fprintf('Zakończono.\n');