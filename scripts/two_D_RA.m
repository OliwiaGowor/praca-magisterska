% =========================================================================
% Title: Chaotic Map 2D-RA (Rastrigin-Ackley) - Bifurcation Diagram
% Description: Generates a classical bifurcation diagram and LE analysis with rigid scale.
% =========================================================================
clear; clc; close all;

%% 1. Parameter Configuration
bias = 1e8;             
beta = 1;               
Euler = exp(1);         
x_init = 0.5;           
y_init = 0.5;           

% Simulation settings
alpha_min = 0;
alpha_max = 100;        
d_alpha = 0.1;          
alpha_values = alpha_min:d_alpha:alpha_max;

n_transient = 1000;     
n_collect = 300;        
total_iter = n_transient + n_collect;

% Container preparation
LE_results = zeros(length(alpha_values), 2);
Bif_X = zeros(1, length(alpha_values) * n_collect);
Bif_A = zeros(1, length(alpha_values) * n_collect);

fprintf('Obliczanie... Trwa generowanie danych dla alpha [0, 100].\n');

%% 2. Main computational loop
idx_bif = 1;
for k = 1:length(alpha_values)
    alpha = alpha_values(k);
    x = x_init; y = y_init;
    Q = eye(2);
    L_sum = zeros(2, 1);
    
    K1 = beta + bias;
    K2 = alpha + bias;
    
    for i = 1:total_iter
        % Jacobian
        r_sq = 0.5 * x^2 + 0.5 * y^2;
        sqrt_term = max(sqrt(r_sq), eps);
        exp1 = exp(-0.2 * sqrt_term);
        exp2 = exp(0.5 * (cos(2*pi*x) + cos(2*pi*y)));
        
        d_exp1_arg_dx = -0.1 * x / sqrt_term;
        d_exp1_arg_dy = -0.1 * y / sqrt_term;
        d_exp2_arg_dx = -pi * sin(2*pi*x);
        d_exp2_arg_dy = -pi * sin(2*pi*y);
        
        J11 = 2*x + 2*pi*K1*sin(2*pi*x) - K2 * exp1 * d_exp1_arg_dx;
        J12 = -K2 * exp1 * d_exp1_arg_dy;
        J21 = -20 * K2 * exp2 * d_exp2_arg_dx;
        J22 = 2*y + 2*pi*K1*sin(2*pi*y) - 20 * K2 * exp2 * d_exp2_arg_dy;
        
        J = [J11, J12; J21, J22];
        [Q, R] = qr(J * Q);
        L_sum = L_sum + log(abs(diag(R)));
        
        % System iteration
        x_next = mod(x^2 - K1*cos(2*pi*x) - K2*exp1, 1);
        y_next = mod(y^2 - K1*cos(2*pi*y) - 20*K2*exp2 + Euler, 1);
        x = x_next; y = y_next;
        
        if i > n_transient
            Bif_X(idx_bif) = x;
            Bif_A(idx_bif) = alpha;
            idx_bif = idx_bif + 1;
        end
    end
    LE_results(k, :) = L_sum / total_iter;
end

%% 3. Plotting Bifurcation Diagram
h_bif = figure('Name', 'Diagram Bifurkacyjny', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
plot(Bif_A, Bif_X, 'k.', 'MarkerSize', 0.5); 
xlabel('\alpha'); ylabel('x');
title('Diagram Bifurkacyjny - Chaotyczne odwzorowanie 2D-RA');
xlim([alpha_min, alpha_max]); ylim([0, 1]);
grid on;

% Saving bifurcation diagram
exportgraphics(h_bif, 'Bifurkacja_2D_RA.png', 'Resolution', 300);
fprintf('Zapisano: Bifurkacja_2D_RA.png\n');

%% 4. Plotting Lyapunov Exponents (Set scale 18-24)
h_le = figure('Name', 'Wykładniki Lapunowa', 'Color', 'w', 'Position', [100, 100, 1000, 450]);
hold on;

plot(alpha_values, LE_results(:, 1), 'b-', 'LineWidth', 1.5, 'DisplayName', '\lambda_1');
plot(alpha_values, LE_results(:, 2), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2, 'DisplayName', '\lambda_2');

ylim([18, 24]); 
xlim([alpha_min, alpha_max]);
xlabel('\alpha'); ylabel('\lambda');
title('Wykładniki Lapunowa - Chaotyczne odwzorowanie 2D-RA');
legend('Location', 'northeast');
grid on; box on;

% Saving Lyapunov exponents
exportgraphics(h_le, 'Lapunov_2D_RA.png', 'Resolution', 300);
fprintf('Zapisano: Lapunov_2D_RA.png\n');

fprintf('Proces zakończony pomyślnie.\n');
