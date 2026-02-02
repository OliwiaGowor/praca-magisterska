% =========================================================================
% Skrypt: Chaotyczne odwzorowanie 2D-RA (Rastrigin-Ackley)
% Opis: Generuje diagram bifurkacyjny i wykres wykładników Lapunowa
%       zgodnie z artykułem arXiv:2509.06754v3.
% =========================================================================

clear; clc; close all;

%% 1. Ustawienia parametrów (zgodnie z sekcją 2.1.2 i Fig. 1, Fig. 3)
% Stałe z artykułu
bias = 1e8;         % Offset (sekcja 2.1.2)
beta = 1;           % Parametr beta stały dla tego eksperymentu
x_init = 0.5;       % Warunek początkowy x
y_init = 0.5;       % Warunek początkowy y
Euler = exp(1);     % Stała e

% Ustawienia symulacji
alpha_values = 0:0.5:100; % Zakres alpha (zgodnie z Fig. 3)
n_iter = 1000;            % Liczba iteracji do odrzucenia (stan przejściowy)
n_keep = 500;             % Liczba punktów do diagramu bifurkacyjnego
total_iter = n_iter + n_keep;

% Przygotowanie tablic na wyniki
LE_history = zeros(length(alpha_values), 2); % Na dwa wykładniki Lapunowa
Bifurcation_x = []; % Na punkty diagramu bifurkacyjnego
Bifurcation_alpha = [];

disp('Rozpoczynanie obliczeń. To może chwilę potrwać...');

%% 2. Główna pętla po parametrze Alpha
for k = 1:length(alpha_values)
    alpha = alpha_values(k);
    
    % Reset zmiennych
    x = x_init;
    y = y_init;
    
    % Inicjalizacja do obliczeń LE (metoda QR)
    Q = eye(2); 
    L_sum = zeros(2, 1);
    
    % Parametry pomocnicze dla czytelności równań
    K1 = beta + bias;
    K2 = alpha + bias;
    
    for i = 1:total_iter
        % --- Obliczanie Jacobianu (do LE) ---
        % Przed aktualizacją x i y musimy policzyć pochodne cząstkowe
        
        % Termy pomocnicze do pochodnych
        sqrt_term = sqrt(0.5*x^2 + 0.5*y^2);
        % Zabezpieczenie przed dzieleniem przez 0
        if sqrt_term == 0, sqrt_term = eps; end
        
        exp_term1 = exp(-0.2 * sqrt_term);
        exp_term2 = exp(0.5 * (cos(2*pi*x) + cos(2*pi*y)));
        
        % Elementy Jacobianu J = [df1/dx, df1/dy; df2/dx, df2/dy]
        % Pochodne f1 (równanie x)
        % d(-0.2 * sqrt(0.5x^2 + 0.5y^2)) / dx = -0.2 * (0.5*2x)/(2*sqrt) = -0.1*x/sqrt
        d_exp1_dx = exp_term1 * (-0.1 * x / sqrt_term);
        d_exp1_dy = exp_term1 * (-0.1 * y / sqrt_term);
        
        J11 = 2*x + 2*pi*K1*sin(2*pi*x) - K2 * d_exp1_dx;
        J12 = -K2 * d_exp1_dy;
        
        % Pochodne f2 (równanie y)
        d_exp2_dx = exp_term2 * 0.5 * (-2*pi*sin(2*pi*x));
        d_exp2_dy = exp_term2 * 0.5 * (-2*pi*sin(2*pi*y));
        
        J21 = -K2 * d_exp2_dx; 
        J22 = 2*y + 2*pi*K1*sin(2*pi*y) - K2 * d_exp2_dy;
        
        J = [J11, J12; J21, J22];
        
        % Aktualizacja macierzy ortogonalnej Q i sumy logarytmów (Algorytm QR)
        [Q, R] = qr(J * Q);
        L_sum = L_sum + log(abs(diag(R)));
        
        % --- Iteracja odwzorowania (Równanie 5 i 6) ---
        x_next = mod(x^2 - K1*cos(2*pi*x) - K2*exp_term1, 1);
        y_next = mod(y^2 - K1*cos(2*pi*y) - K2*exp_term2 + Euler, 1);
        
        x = x_next;
        y = y_next;
        
        % --- Zapisywanie danych do diagramu bifurkacyjnego ---
        if i > n_iter
            Bifurcation_x(end+1) = x; %#ok<SAGROW>
            Bifurcation_alpha(end+1) = alpha; %#ok<SAGROW>
        end
    end
    
    % Średnia dla LE
    LE_history(k, :) = L_sum / total_iter;
end

disp('Obliczenia zakończone. Generowanie wykresów...');

%% 3. Rysowanie i zapisywanie Diagramu Bifurkacyjnego
figure('Name', 'Diagram Bifurkacyjny 2D-RA', 'Color', 'w');
% Używamy małych kropek ('.') i lekkiej przezroczystości (jeśli wersja wspiera),
% aby symulować efekt "heatmapy" gęstości punktów jak na Fig. 1.
scatter(Bifurcation_alpha, Bifurcation_x, 1, 'b', '.');
xlim([min(alpha_values), max(alpha_values)]);
ylim([0, 1]);
xlabel('\alpha');
ylabel('x');
title('Diagram Bifurkacyjny funkcji 2D-RA');
box on;
grid on;

% Zapis do pliku
print('Chaotic_2D_RA_Bifurcation', '-dpng', '-r300');
disp('Zapisano: Chaotic_2D_RA_Bifurcation.png');

%% 4. Rysowanie i zapisywanie Wykładników Lapunowa (LE)
figure('Name', 'Wykładniki Lapunowa 2D-RA', 'Color', 'w');
plot(alpha_values, LE_history(:, 1), 'b-', 'LineWidth', 1.5); hold on;
plot(alpha_values, LE_history(:, 2), 'r--', 'LineWidth', 1.5);
xlim([min(alpha_values), max(alpha_values)]);
% Opcjonalnie dostosuj ylim do zakresu z artykułu (ok. 19-22)
xlabel('\alpha');
ylabel('Wykładniki Lapunowa (LE)');
legend('\lambda_1', '\lambda_2');
title('Analiza Wykładników Lapunowa funkcji 2D-RA');
box on;
grid on;

% Zapis do pliku
print('Chaotic_2D_RA_Lyapunov', '-dpng', '-r300');
disp('Zapisano: Chaotic_2D_RA_Lyapunov.png');

disp('Gotowe.');