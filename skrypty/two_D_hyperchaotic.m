%% Parametry i inicjalizacja systemu (zgodnie z sekcją 2 artykułu [cite: 59, 64])
clear; clc;

% Nazwa układu do nazewnictwa plików i tytułów
sys_name = 'Nowy_dwuwymiarowy_uklad_hiperchaotyczny';
sys_title = 'Nowy dwuwymiarowy układ hiperchaotyczny';

% Parametry kontrolne i zakresy
a_range = linspace(0.1, 60, 1000); % Zakres parametru 'a' [cite: 81]
b = 30;                            % Stała wartość 'b' dla diagramu [cite: 81]
x0 = 0.2; y0 = 0.3;                % Wartości początkowe 
n_transient = 500;                 % Iteracje przejściowe do pominięcia
n_iter = 1000;                     % Liczba iteracji do analizy statystycznej

% Inicjalizacja struktur danych
num_bins = 500;                    % Rozdzielczość pionowa heatmapy
heatmap_matrix = zeros(num_bins, length(a_range));
x_edges = linspace(0, 1, num_bins + 1);
LE1 = zeros(size(a_range));
LE2 = zeros(size(a_range));

%% Pętla obliczeniowa - Dynamika systemu i Wykładniki Lapunowa
% Obliczenia oparte na równaniach: 
% x(n+1) = sin^2(a*pi*x(n) + b*y(n)) [cite: 67]
% y(n+1) = cos^2(b*pi/y(n) + a*x(n)) [cite: 68]

for i = 1:length(a_range)
    a = a_range(i);
    x = x0; y = y0;
    
    % 1. Faza przejściowa (stabilizacja trajektorii)
    for n = 1:n_transient
        x_next = sin(a*pi*x + b*y)^2; 
        y_next = cos(b*pi/y + a*x)^2; 
        x = x_next; y = y_next;
    end
    
    % 2. Zbieranie danych do heatmapy i obliczanie LE
    current_x_series = zeros(1, n_iter);
    sum_log_eig = zeros(2, 1);
    
    for n = 1:n_iter
        % Obliczenie kolejnego kroku systemu [cite: 67, 68]
        x_next = sin(a*pi*x + b*y)^2;
        y_next = cos(b*pi/y + a*x)^2;
        
        % Macierz Jacobiego do obliczenia wykładników Lapunowa [cite: 164]
        U = a*pi*x + b*y;
        V = b*pi/y + a*x;
        
        % Elementy macierzy Jacobiego (pochodne cząstkowe)
        J11 = sin(2*U) * a * pi;
        J12 = sin(2*U) * b;
        J21 = -sin(2*V) * a;
        J22 = -sin(2*V) * (-b*pi / (y^2));
        
        % Obliczanie wartości własnych dla oszacowania chaosu [cite: 172]
        eig_vals = abs(eig([J11, J12; J21, J22]));
        sum_log_eig = sum_log_eig + log(max(eig_vals, 1e-10));
        
        x = x_next; y = y_next;
        current_x_series(n) = x;
    end
    
    % Zapis histogramu gęstości dla bieżącego 'a'
    heatmap_matrix(:, i) = histcounts(current_x_series, x_edges)';
    
    % Obliczanie średnich wykładników Lapunowa [cite: 172]
    LE_vals = sort(sum_log_eig / n_iter, 'descend');
    LE1(i) = LE_vals(1);
    LE2(i) = LE_vals(2);
end

%% Generowanie i zapisywanie wykresów do plików

% --- Wykres 1: Heatmapa Bifurkacyjna ---
fig1 = figure('Color', 'w', 'Name', 'Diagram Bifurkacyjny');
imagesc(a_range, [0 1], log1p(heatmap_matrix)); 
set(gca, 'YDir', 'normal');
colormap(jet);
cb = colorbar;
cb.Label.String = 'Wykładniki Lapunowa';
title(['Diagram bifurkacyjny - ', sys_title]);
xlabel('Parametr kontrolny a');
ylabel('Zmienna stanu x_n');
grid on;

% Zapis pliku z nową nazwą
file_name1 = sprintf('%s_bifurkacja.png', sys_name);
saveas(fig1, file_name1);

% --- Wykres 2: Wykładniki Lapunowa ---
fig2 = figure('Color', 'w', 'Name', 'Wykladniki Lapunowa');
plot(a_range, LE1, 'r', 'LineWidth', 1.2); hold on;
plot(a_range, LE2, 'b', 'LineWidth', 1.2);
line([0 60], [0 0], 'Color', 'k', 'LineStyle', '--'); % Linia progu chaosu
title(['Wykładniki Lapunowa - ', sys_title]);
xlabel('Parametr kontrolny a');
ylabel('Wartości LE');
legend('LE_1', 'LE_2', 'Location', 'best');
grid on;
xlim([0 60]);

% Zapis pliku z nową nazwą
file_name2 = sprintf('%s_le.png', sys_name);
saveas(fig2, file_name2);

fprintf('Zapisano pliki:\n1. %s\n2. %s\n', file_name1, file_name2);