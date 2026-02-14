% Two-stage Logistic Mapping
% Code generates bifurcation diagram and Lyapunov exponent plot.

clear; clc; close all;

%% 1. Simulation parameters
gamma_values = 2.9:0.001:4.0;

% Number of iterations
n_total = 2000;       % Total number of steps
n_transient = 1000;   % Discarded initial steps (transient state)
n_plot = n_total - n_transient; % Number of points for the plot

% Initial value
x0 = 0.41;

% Data vector preparation
bif_gamma = [];
bif_x = [];
lyap_exp = zeros(size(gamma_values));

fprintf('Obliczanie diagramu bifurkacyjnego i wykładnika Lapunowa...\n');

%% 2. Computational loop
for i = 1:length(gamma_values)
    g = gamma_values(i);
    x = x0;
    
    L_sum = 0; % Sum for Lyapunov exponent
    
    % Initial iteration (to establish attractor)
    for k = 1:n_transient
        if x < 0.5
            % Branch 1: Left parabola
            x = 4 * g * x * (0.5 - x);
        else
            x = 4 * g * (x - 0.5) * (1 - x);
        end
    end
    
    % Collecting data for plots
    x_orbit = zeros(1, n_plot);
    for k = 1:n_plot
        % Calculation of value and derivative
        if x < 0.5
            x_next = 4 * g * x * (0.5 - x);
            df = 4 * g * (0.5 - 2 * x); % Derivative
        else
            x_next = 4 * g * (x - 0.5) * (1 - x);
            % Derivative for 4*g*(-x^2 + 1.5x - 0.5)
            df = 4 * g * (1.5 - 2 * x);
        end
        x = x_next;
        x_orbit(k) = x;
        
        % Summing derivative logarithm (formula 2.1 from the article)
        L_sum = L_sum + log(abs(df));
    end
    
    % Saving results
    bif_gamma = [bif_gamma, g * ones(size(x_orbit))];
    bif_x = [bif_x, x_orbit];
    lyap_exp(i) = L_sum / n_plot;
end

%% 3. Plotting and saving Bifurcation Diagram
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

%% 4. Plotting and saving Lyapunov Exponent
figure('Name', 'Wykładnik Lapunowa', 'Color', 'w', 'Visible', 'off');
plot(gamma_values, lyap_exp, 'b-', 'LineWidth', 1.2);
hold on;
plot([2.9, 4.0], [0, 0], 'k--'); % Zero line
xlim([2.9, 4.0]);
% Automatic Y-axis scaling
ylim([min(lyap_exp)-0.5, max(lyap_exp)+0.5]);
xlabel('\gamma');
ylabel('\lambda');
title('Wykładnik Lapunowa - Dwuetapowe odwzorowanie logistyczne');
grid on;
saveas(gcf, 'dwuetapowe_odwzorowanie_logistyczne_lyapunov.png');
fprintf('Zapisano: dwuetapowe_odwzorowanie_logistyczne_lyapunov.png\n');

fprintf('Zakończono pomyślnie.\n');
