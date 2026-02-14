% =========================================================================
% Script Name: One-dimensional Cosine-Polynomial Map
% Description: Generates bifurcation diagram and Lyapunov exponent plot 
% for the chaotic map from the article.
% =========================================================================

clear; clc; close all;

%% 1. Simulation Parameters
mu_min = 0;
mu_max = 5;
n_points = 1000;            % Parameter mu resolution
mu_values = linspace(mu_min, mu_max, n_points);

a = 4;                      % Constant value of parameter a
x0 = 0.1;                   % Initial condition

iterations = 1000;          % Total iterations per point
transient = 500;            % Number of transient iterations (discarded)

% Data vector preparation
bifurcation_x = [];
bifurcation_mu = [];
lyapunov_exponents = zeros(1, n_points);

%% 2. Main Computational Loop
fprintf('Trwa generowanie danych dla zakresu mu = [%.1f, %.1f]...\n', mu_min, mu_max);

for k = 1:n_points
    mu = mu_values(k);
    x = x0;
    
    % Variable for summing logarithm of derivatives for LE
    log_sum = 0;
    
    % Transient iterations
    for i = 1:transient
        % Map equation (2): x_k = cos(mu * (x^3 + x) + a) 
        x = cos(mu * (x^3 + x) + a);
    end
    
    % Actual iterations for diagram and LE
    current_orbit = zeros(1, iterations - transient);
    
    for i = 1:(iterations - transient)
        x_prev = x;
        
        % 1. Calculation of new x value
        x = cos(mu * (x_prev^3 + x_prev) + a);
        
        % Saving orbit point
        current_orbit(i) = x;
        
        % 2. Derivative calculation for Lyapunov Exponent
        % f(x) = cos(u), where u = mu*(x^3 + x) + a
        % f'(x) = -sin(u) * u'
        % u' = mu * (3*x^2 + 1)
        derivative = -sin(mu * (x_prev^3 + x_prev) + a) * mu * (3 * x_prev^2 + 1);
        
        % Summing logarithm of derivative magnitude
        log_sum = log_sum + log(abs(derivative));
    end
    
    % Saving data for bifurcation diagram
    bifurcation_x = [bifurcation_x, current_orbit];
    bifurcation_mu = [bifurcation_mu, repmat(mu, 1, length(current_orbit))];
    
    % Calculating average LE (standard formula for 1D maps)
    lyapunov_exponents(k) = log_sum / (iterations - transient);
end

%% 3. Plotting

% 1. Bifurcation Diagram
fig1 = figure('Name', 'Diagram Bifurkacyjny', 'Color', 'w');
% Using dots with high transparency (if Matlab version supports) 
% or small size to simulate point density "heatmap".
plot(bifurcation_mu, bifurcation_x, '.', 'MarkerSize', 1, 'Color', [0 0 0]);
xlim([mu_min mu_max]);
ylim([-1 1]);
xlabel('\mu');
ylabel('x');
title(['Diagram Bifurkacyjny (a = ' num2str(a) ')']);
grid on;

% 2. Lyapunov Exponent
fig2 = figure('Name', 'Wykładnik Lapunowa', 'Color', 'w');
plot(mu_values, lyapunov_exponents, 'b-', 'LineWidth', 1.2);
hold on;
plot([mu_min, mu_max], [0, 0], 'r--'); % Zero reference line
xlim([mu_min mu_max]);
ylim([min(lyapunov_exponents)-0.5, max(lyapunov_exponents)+0.5]);
xlabel('\mu');
ylabel('Wykładnik Lapunowa (LE)');
title(['Wykładnik Lapunowa (a = ' num2str(a) ')']);
legend('LE', 'Granica chaosu (0)');
grid on;

%% 4. Saving to files
filename_base = 'Jednowymiarowe_odwzorowanie_kosinusowo_wielomianowe';

% Saving bifurcation diagram
file_bif = [filename_base '_Bifurkacja.png'];
saveas(fig1, file_bif);
fprintf('Zapisano diagram bifurkacyjny do: %s\n', file_bif);

% Saving Lyapunov plot
file_lap = [filename_base '_Lapunov.png'];
saveas(fig2, file_lap);
fprintf('Zapisano wykres Lapunowa do: %s\n', file_lap);

fprintf('Gotowe.\n');
