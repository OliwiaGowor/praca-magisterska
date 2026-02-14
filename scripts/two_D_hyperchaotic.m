% =========================================================================
% Title: New 2D Hyperchaotic System Analysis
% Description: Generates bifurcation diagram and Lyapunov exponents.
% =========================================================================
clear; clc; close all;

%% 1. Parameter Configuration
a_range = linspace(0, 60, 600); 
b = 30;                         
x0 = 0.2;                        
y0 = 0.3;
n_transient = 200;                
n_iter = 400;                     

LE1 = zeros(size(a_range));
LE2 = zeros(size(a_range));
sys_title = 'Nowy dwuwymiarowy układ hiperchaotyczny';

%% 2. Generating and Saving Bifurcation Diagram
fig1 = figure('Name', 'Bifurkacja', 'Position', [100, 100, 800, 500]);
hold on;

disp('Generowanie diagramu bifurkacyjnego...');
for i = 1:length(a_range)
    a = a_range(i);
    x = x0; y = y0;
    
    % Transient iterations
    for n = 1:n_transient
        x_next = sin(a*pi*x + b*y)^2;
        y_next = cos(b*pi/y + a*x)^2;
        x = x_next; y = y_next;
    end
    
    % Data collection
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

% Saving bifurcation diagram
saveas(fig1, 'diagram_bifurkacyjny.png');
disp('Zapisano: diagram_bifurkacyjny.png');

%% 3. Calculating and Saving Lyapunov Exponents
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
        
        % Jacobian matrix calculation
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

% Saving Lyapunov plot
saveas(fig2, 'wykladniki_lapunowa.png');
disp('Zapisano: wykladniki_lapunowa.png');

%% 4. Saving Numerical Data

% Preparing LE results table
resultsTable = table(a_range', LE1', LE2', ...
    'VariableNames', {'Parametr_a', 'LE1', 'LE2'});

% Saving data to files (added to match the print statement)
writetable(resultsTable, 'dane_lapunowa.csv');
save('pelne_dane_systemu.mat', 'a_range', 'LE1', 'LE2', 'b', 'x0', 'y0');

fprintf('\nZakończono sukcesem!\nUtworzono pliki graficzne:\n - diagram_bifurkacyjny.png\n - wykladniki_lapunowa.png\nOraz pliki danych:\n - dane_lapunowa.csv\n - pelne_dane_systemu.mat\n');
