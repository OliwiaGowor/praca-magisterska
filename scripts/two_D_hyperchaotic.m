%% Parametry i inicjalizacja
clear; clc; close all; % close all zamyka poprzednie okna
a_range = linspace(0, 60, 600); 
b = 30;                          
x0 = 0.2;                        
y0 = 0.3;
n_transient = 200;               
n_iter = 400;                    

LE1 = zeros(size(a_range));
LE2 = zeros(size(a_range));
sys_title = 'Nowy dwuwymiarowy układ hiperchaotyczny';

%% 1. Generowanie i zapisywanie Diagramu Bifurkacyjnego
% Tworzymy osobne okno dla pierwszego wykresu
fig1 = figure('Name', 'Bifurkacja', 'Position', [100, 100, 800, 500]);
hold on;

disp('Generowanie diagramu bifurkacyjnego...');
for i = 1:length(a_range)
    a = a_range(i);
    x = x0; y = y0;
    % Iteracje przejściowe
    for n = 1:n_transient
        x_next = sin(a*pi*x + b*y)^2;
        y_next = cos(b*pi/y + a*x)^2;
        x = x_next; y = y_next;
    end
    % Zbieranie danych
    x_data = zeros(1, n_iter);
    for n = 1:n_iter
        x_next = sin(a*pi*x + b*y)^2;
        y_next = cos(b*pi/y + a*x)^2;
        x = x_next; y = y_next;
        x_data(n) = x;
    end
    plot(ones(1, n_iter)*a, x_data, 'k.', 'MarkerSize', 0.5);
end

title(['Diagram bifurkacyjny - ', sys_title]);
xlabel('a'); ylabel('x');
xlim([0 60]); ylim([0 1]);
grid on;

% ZAPIS DIAGRAMU BIFURKACYJNEGO
saveas(fig1, 'diagram_bifurkacyjny.png');
disp('Zapisano: diagram_bifurkacyjny.png');


%% 2. Obliczanie i zapisywanie Wykładników Lapunowa
% Tworzymy osobne okno dla drugiego wykresu
fig2 = figure('Name', 'Wykladniki Lapunowa', 'Position', [150, 150, 800, 500]);

disp('Obliczanie wykładników Lapunowa...');
N_LE = 1000; 
for i = 1:length(a_range)
    a = a_range(i);
    x = x0; y = y0;
    sum_log_eig = zeros(2, 1);
    
    for n = 1:N_LE
        x_next = sin(a*pi*x + b*y)^2;
        y_next = cos(b*pi/y + a*x)^2;
        
        J11 = 2*sin(a*pi*x + b*y)*cos(a*pi*x + b*y)*(a*pi);
        J12 = 2*sin(a*pi*x + b*y)*cos(a*pi*x + b*y)*(b);
        J21 = 2*cos(b*pi/y + a*x)*(-sin(b*pi/y + a*x))*(a);
        J22 = 2*cos(b*pi/y + a*x)*(-sin(b*pi/y + a*x))*(-b*pi/(y^2));
        
        Jac = [J11, J12; J21, J22];
        x = x_next; y = y_next;
        
        eig_vals = abs(eig(Jac));
        sum_log_eig = sum_log_eig + log(max(eig_vals, 1e-10)); 
    end
    
    LE_vals = sort(sum_log_eig / N_LE, 'descend');
    LE1(i) = LE_vals(1);
    LE2(i) = LE_vals(2);
end

plot(a_range, LE1, 'r', 'LineWidth', 1); hold on;
plot(a_range, LE2, 'b', 'LineWidth', 1);
title(['Wykładniki Lapunowa - ', sys_title]);
xlabel('a'); ylabel('\lambda');
legend('LE_1', 'LE_2', 'Location', 'southeast');
grid on;
xlim([0 60]);

% ZAPIS WYKRESU LAPUNOWA
saveas(fig2, 'wykladniki_lapunowa.png');
disp('Zapisano: wykladniki_lapunowa.png');


%% 3. ZAPISYWANIE DANYCH LICZBOWYCH

% Przygotowanie tabeli z wynikami LE
resultsTable = table(a_range', LE1', LE2', ...
    'VariableNames', {'Parametr_a', 'LE1', 'LE2'});

% Zapis do pliku CSV
writetable(resultsTable, 'dane_lapunowa.csv');

% Zapis wszystkich zmiennych do .mat
save('pelne_dane_systemu.mat');

fprintf('\nZakończono sukcesem!\nUtworzono pliki graficzne:\n - diagram_bifurkacyjny.png\n - wykladniki_lapunowa.png\nOraz pliki danych:\n - dane_lapunowa.csv\n - pelne_dane_systemu.mat\n');