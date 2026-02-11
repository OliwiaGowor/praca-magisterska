% Dwuetapowe Odwzorowanie Logistyczne (Two-stage Logistic Mapping)
% Kod generuje diagram bifurkacyjny i wykres wykładnika Lapunowa.
% Zastosowano korektę wzoru z artykułu, aby uniknąć błędów ucieczki do nieskończoności.

clear; clc; close all;

%% 1. Parametry symulacji
% Zakres parametru gamma (zgodnie z Rysunkiem 2 w artykule)
gamma_values = 2.9:0.001:4.0;

% Liczba iteracji
n_total = 2000;       % Całkowita liczba kroków
n_transient = 1000;   % Odrzucane kroki początkowe (stan przejściowy)
n_plot = n_total - n_transient; % Liczba punktów do wykresu

% Wartość początkowa (z sekcji eksperymentalnej artykułu)
x0 = 0.41;

% Przygotowanie wektorów danych
bif_gamma = [];
bif_x = [];
lyap_exp = zeros(size(gamma_values));

fprintf('Obliczanie diagramu bifurkacyjnego i wykładnika Lapunowa...\n');

%% 2. Pętla obliczeniowa
for i = 1:length(gamma_values)
    g = gamma_values(i);
    x = x0;
    
    L_sum = 0; % Suma do wykładnika Lapunowa
    
    % Iteracja wstępna (do ustalenia atraktora)
    for k = 1:n_transient
        % Zastosowanie poprawionego wzoru symetrycznego
        if x < 0.5
            % Gałąź 1: Lewa parabola
            x = 4 * g * x * (0.5 - x);
        else
            % Gałąź 2: Prawa parabola (Symetryczna)
            % UWAGA: Oryginalny wzór z artykułu [1 - 4*g*x*(x-0.5)] jest błędny
            % i powoduje ucieczkę do nieskończoności dla g > 0.5.
            % Poniższy wzór to poprawna forma 'Two-stage':
            x = 4 * g * (x - 0.5) * (1 - x);
        end
    end
    
    % Zbieranie danych do wykresów
    x_orbit = zeros(1, n_plot);
    for k = 1:n_plot
        % Obliczenie wartości i pochodnej
        if x < 0.5
            x_next = 4 * g * x * (0.5 - x);
            df = 4 * g * (0.5 - 2 * x); % Pochodna
        else
            x_next = 4 * g * (x - 0.5) * (1 - x);
            % Pochodna dla 4*g*(-x^2 + 1.5x - 0.5)
            df = 4 * g * (1.5 - 2 * x);
        end
        x = x_next;
        x_orbit(k) = x;
        
        % Sumowanie logarytmu pochodnej (wzór 2.1 z artykułu)
        L_sum = L_sum + log(abs(df));
    end
    
    % Zapis wyników
    bif_gamma = [bif_gamma, g * ones(size(x_orbit))];
    bif_x = [bif_x, x_orbit];
    lyap_exp(i) = L_sum / n_plot;
end

%% 3. Rysowanie i zapisywanie Diagramu Bifurkacyjnego
figure('Name', 'Diagram Bifurkacyjny', 'Color', 'w', 'Visible', 'off');
plot(bif_gamma, bif_x, '.', 'MarkerSize', 1, 'Color', 'k');
xlim([2.9, 4.0]);
ylim([0, 1]);
xlabel('\gamma');
ylabel('x');
title('Diagram Bifurkacyjny - Dwuetapowe odwzorowanie logistyczne');
grid on;
saveas(gcf, 'dwuetapowe_odwzorowanie_logistyczne_bifurkacja.png');
fprintf('Zapisano: dwuetapowe_odwzorowanie_logistyczne_bifurkacja.png\n');

%% 4. Rysowanie i zapisywanie Wykładnika Lapunowa
figure('Name', 'Wykładnik Lapunowa', 'Color', 'w', 'Visible', 'off');
plot(gamma_values, lyap_exp, 'b-', 'LineWidth', 1.2);
hold on;
plot([2.9, 4.0], [0, 0], 'k--'); % Linia zerowa
xlim([2.9, 4.0]);
% Automatyczne skalowanie osi Y
ylim([min(lyap_exp)-0.5, max(lyap_exp)+0.5]);
xlabel('\gamma');
ylabel('\lambda');
title('Wykładnik Lapunowa - Dwuetapowe odwzorowanie logistyczne');
grid on;
saveas(gcf, 'dwuetapowe_odwzorowanie_logistyczne_lyapunov.png');
fprintf('Zapisano: dwuetapowe_odwzorowanie_logistyczne_lyapunov.png\n');

fprintf('Zakończono pomyślnie.\n');