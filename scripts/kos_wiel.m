% Nazwa skryptu: Jednowymiarowe_odwzorowanie_kosinusowo_wielomianowe.m
% Opis: Generuje diagram bifurkacyjny i wykres wykładnika Lapunowa 
% dla mapy chaotycznej z artykułu.

clear; clc; close all;

%% Parametry symulacji
% Zgodnie z sekcją 2 artykułu (Fig. 2 i Fig. 5) [cite: 165, 232]
mu_min = 0;
mu_max = 5;
n_points = 1000;            % Rozdzielczość parametru mu
mu_values = linspace(mu_min, mu_max, n_points);

a = 4;                      % Stała wartość parametru a [cite: 133]
x0 = 0.1;                   % Warunek początkowy [cite: 133]

iterations = 1000;          % Całkowita liczba iteracji na punkt
transient = 500;            % Liczba iteracji przejściowych (odrzucanych)

% Przygotowanie wektorów na dane
bifurcation_x = [];
bifurcation_mu = [];
lyapunov_exponents = zeros(1, n_points);

%% Główna pętla obliczeniowa
fprintf('Trwa generowanie danych dla zakresu mu = [%.1f, %.1f]...\n', mu_min, mu_max);

for k = 1:n_points
    mu = mu_values(k);
    x = x0;
    
    % Zmienna do sumowania logarytmów pochodnych dla LE
    log_sum = 0;
    
    % Iteracje przejściowe (transient)
    for i = 1:transient
        % Równanie mapy (2): x_k = cos(mu * (x^3 + x) + a) 
        x = cos(mu * (x^3 + x) + a);
    end
    
    % Właściwe iteracje do diagramu i LE
    current_orbit = zeros(1, iterations - transient);
    
    for i = 1:(iterations - transient)
        x_prev = x;
        
        % 1. Obliczenie nowej wartości x
        x = cos(mu * (x_prev^3 + x_prev) + a);
        
        % Zapisanie punktu orbity
        current_orbit(i) = x;
        
        % 2. Obliczenie pochodnej do Wykładnika Lapunowa
        % f(x) = cos(u), gdzie u = mu*(x^3 + x) + a
        % f'(x) = -sin(u) * u'
        % u' = mu * (3*x^2 + 1)
        derivative = -sin(mu * (x_prev^3 + x_prev) + a) * mu * (3 * x_prev^2 + 1);
        
        % Sumowanie logarytmu modułu pochodnej
        log_sum = log_sum + log(abs(derivative));
    end
    
    % Zapisanie danych do diagramu bifurkacyjnego
    bifurcation_x = [bifurcation_x, current_orbit];
    bifurcation_mu = [bifurcation_mu, repmat(mu, 1, length(current_orbit))];
    
    % Obliczenie średniego LE (wzór standardowy dla map 1D)
    lyapunov_exponents(k) = log_sum / (iterations - transient);
end

%% Rysowanie Wykresów

% 1. Diagram Bifurkacyjny
fig1 = figure('Name', 'Diagram Bifurkacyjny', 'Color', 'w');
% Używamy kropek o dużej przezroczystości (jeśli wersja Matlaba obsługuje) 
% lub małego rozmiaru, aby symulować "heatmapę" gęstości punktów.
plot(bifurcation_mu, bifurcation_x, '.', 'MarkerSize', 1, 'Color', [0 0 0]);
xlim([mu_min mu_max]);
ylim([-1 1]);
xlabel('\mu');
ylabel('x');
title(['Diagram Bifurkacyjny (a = ' num2str(a) ')']);
grid on;

% 2. Wykładnik Lapunowa
fig2 = figure('Name', 'Wykładnik Lapunowa', 'Color', 'w');
plot(mu_values, lyapunov_exponents, 'b-', 'LineWidth', 1.2);
hold on;
plot([mu_min, mu_max], [0, 0], 'r--'); % Linia odniesienia zero
xlim([mu_min mu_max]);
ylim([min(lyapunov_exponents)-0.5, max(lyapunov_exponents)+0.5]);
xlabel('\mu');
ylabel('Wykładnik Lapunowa (LE)');
title(['Wykładnik Lapunowa (a = ' num2str(a) ')']);
legend('LE', 'Granica chaosu (0)');
grid on;

%% Zapisywanie do plików
filename_base = 'Jednowymiarowe_odwzorowanie_kosinusowo_wielomianowe';

% Zapis diagramu bifurkacyjnego
file_bif = [filename_base '_Bifurkacja.png'];
saveas(fig1, file_bif);
fprintf('Zapisano diagram bifurkacyjny do: %s\n', file_bif);

% Zapis wykresu Lapunowa
file_lap = [filename_base '_Lapunov.png'];
saveas(fig2, file_lap);
fprintf('Zapisano wykres Lapunowa do: %s\n', file_lap);

fprintf('Gotowe.\n');