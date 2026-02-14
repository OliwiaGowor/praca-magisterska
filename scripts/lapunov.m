clear; clc; close all;

%% === GLOBAL SETTINGS ===
res_x = 800;          % Horizontal resolution (number of parameter steps)
n_discard = 500;      % Discarding transient state
n_keep = 400;         % Number of points to plot on the diagram
n_iter_lle = 1000;    % Number of iterations to calculate LLE

if isempty(gcp('nocreate')), parpool; end
fprintf('=== START SYMULACJI (DIAGRAMY CZARNO-BIAŁE) ===\n');

%% --- 1. COSINE-POLYNOMIAL MAP ---
fprintf('[1/7] Układ Kosinusowo-Wielomianowy...\n');
a_cos = 4;
mu_vals = linspace(0, 5, res_x);
bifur_1 = zeros(n_keep, res_x); 
lle_1 = zeros(1, res_x);

parfor i = 1:res_x
    mu = mu_vals(i);
    x = 0.1;
    for k = 1:n_discard
        arg = mu*(x^3+x) + a_cos; x = cos(arg);
    end
    for k = 1:n_keep
        arg = mu*(x^3+x) + a_cos; x = cos(arg);
        bifur_1(k, i) = x;
    end
    
    x=0.1; sum_L=0;
    for k = 1:n_iter_lle+n_discard
        arg = mu*(x^3+x) + a_cos; x_next = cos(arg);
        if k>n_discard
            deriv = -sin(arg)*mu*(3*x^2+1);
            sum_L = sum_L + log(max(abs(deriv), 1e-10));
        end
        x = x_next;
    end
    lle_1(i) = sum_L/n_iter_lle;
end
save_plots(1, mu_vals, bifur_1, lle_1, '1_CosPoly', '\mu', 'Odwzorowanie kosinusowo-wielomianowe', [-1.1 1.1]);


%% --- 2. 3D LOGISTIC MAP ---
fprintf('[2/7] Układ 3D Logistic Map...\n');
beta = 0.017; gamma = 0.0012;
alpha_vals = linspace(3.65, 4.0, res_x);
bifur_2 = zeros(n_keep, res_x); 
lle_2 = zeros(1, res_x);

parfor i = 1:res_x
    alpha = alpha_vals(i);
    x=0.2487; y=0.35; z=0.10001;
    for k=1:n_discard
        xn = alpha*x*(1-x) + beta*y^2*x + gamma*z^3;
        yn = alpha*y*(1-y) + beta*z^2*y + gamma*x^3;
        zn = alpha*z*(1-z) + beta*x^2*z + gamma*y^3;
        x=xn; y=yn; z=zn;
    end
    for k=1:n_keep
        xn = alpha*x*(1-x) + beta*y^2*x + gamma*z^3;
        yn = alpha*y*(1-y) + beta*z^2*y + gamma*x^3;
        zn = alpha*z*(1-z) + beta*x^2*z + gamma*y^3;
        bifur_2(k, i) = xn;
        x=xn; y=yn; z=zn;
    end
    
    x=0.2487; y=0.35; z=0.10001; Q=eye(3); sum_L=0;
    for k=1:n_iter_lle+n_discard
        xn = alpha*x*(1-x) + beta*y^2*x + gamma*z^3;
        yn = alpha*y*(1-y) + beta*z^2*y + gamma*x^3;
        zn = alpha*z*(1-z) + beta*x^2*z + gamma*y^3;
        if k>n_discard
            J = [alpha*(1-2*x)+beta*y^2, 2*beta*x*y, 3*gamma*z^2;
                 3*gamma*x^2, alpha*(1-2*y)+beta*z^2, 2*beta*y*z;
                 2*beta*x*z, 3*gamma*y^2, alpha*(1-2*z)+beta*x^2];
            [Q, R] = qr(J*Q); sum_L = sum_L + log(abs(R(1,1)));
        end
        x=xn; y=yn; z=zn;
    end
    lle_2(i) = sum_L/n_iter_lle;
end
save_plots(2, alpha_vals, bifur_2, lle_2, '2_3DLogistic', '\alpha', 'Trójwymiarowe logistyczne odwzorowanie chaotyczne', [0 1]);


%% --- 3. SOBOLEVA-MODULO MAP ---
fprintf('[3/7] Układ Soboleva-Modulo...\n');
a_sob=5; b_sob=5; K_sob=1; B_sob=1; C_sob=1; D_sob=1;
A_vals = linspace(0, 5, res_x);
bifur_3 = zeros(n_keep, res_x); 
lle_3 = zeros(1, res_x);

parfor i = 1:res_x
    A_val = A_vals(i);
    x=0.1;
    for k=1:n_discard
        u = x/K_sob;
        smht = (exp(A_val*u)-exp(-B_sob*u))/(exp(C_sob*u)+exp(-D_sob*u));
        x = mod(x + a_sob + b_sob*smht, K_sob);
    end
    for k=1:n_keep
        u = x/K_sob;
        smht = (exp(A_val*u)-exp(-B_sob*u))/(exp(C_sob*u)+exp(-D_sob*u));
        x = mod(x + a_sob + b_sob*smht, K_sob);
        bifur_3(k, i) = x;
    end
    
    x=0.1; sum_L=0;
    for k=1:n_iter_lle+n_discard
        u = x/K_sob;
        num = exp(A_val*u)-exp(-B_sob*u); den = exp(C_sob*u)+exp(-D_sob*u);
        smht = num/den;
        xn = mod(x + a_sob + b_sob*smht, K_sob);
        if k>n_discard
             d_num = A_val*exp(A_val*u) + B_sob*exp(-B_sob*u);
             d_den = C_sob*exp(C_sob*u) - D_sob*exp(-D_sob*u);
             smht_p = (d_num*den - num*d_den)/(den^2);
             deriv = 1 + (b_sob/K_sob)*smht_p;
             sum_L = sum_L + log(max(abs(deriv), 1e-10));
        end
        x=xn;
    end
    lle_3(i) = sum_L/n_iter_lle;
end
save_plots(3, A_vals, bifur_3, lle_3, '3_Soboleva', 'A', 'Odwzorowanie Soboleva-Modulo', [0 1]);


%% --- 4. 2D-RA MAP ---
fprintf('[4/7] Układ 2D-RA Map...\n');
bias = 1e8; beta_ra = 1;
alpha_vals_int = 0:(res_x-1); 
bifur_4 = zeros(n_keep, res_x); 
lle_4 = zeros(2, res_x);

parfor i = 1:res_x
    alpha = alpha_vals_int(i);
    const_beta = beta_ra + bias; const_alpha = alpha + bias;
    x=0.5; y=0.5;
    
    % Transient
    for k=1:n_discard
        ts = sqrt(0.5*x^2+0.5*y^2); if ts<1e-10, ts=1e-10; end
        ex = exp(-0.2*ts); ey = exp(0.5*cos(2*pi*x)+0.5*cos(2*pi*y));
        xn = mod(x^2 - const_beta*cos(2*pi*x) - const_alpha*ex, 1);
        yn = mod(y^2 - const_beta*cos(2*pi*y) - const_alpha*ey + exp(1), 1);
        x=xn; y=yn;
    end
    % Recording
    for k=1:n_keep
        ts = sqrt(0.5*x^2+0.5*y^2); if ts<1e-10, ts=1e-10; end
        ex = exp(-0.2*ts); ey = exp(0.5*cos(2*pi*x)+0.5*cos(2*pi*y));
        xn = mod(x^2 - const_beta*cos(2*pi*x) - const_alpha*ex, 1);
        yn = mod(y^2 - const_beta*cos(2*pi*y) - const_alpha*ey + exp(1), 1);
        bifur_4(k, i) = xn;
        x=xn; y=yn;
    end
    
    % LLE
    x=0.5; y=0.5; Q=eye(2); s1=0; s2=0;
    for k=1:n_iter_lle+n_discard
        ts = sqrt(0.5*x^2+0.5*y^2); if ts<1e-10, ts=1e-10; end
        ex = exp(-0.2*ts); ey = exp(0.5*cos(2*pi*x)+0.5*cos(2*pi*y));
        xn = mod(x^2 - const_beta*cos(2*pi*x) - const_alpha*ex, 1);
        yn = mod(y^2 - const_beta*cos(2*pi*y) - const_alpha*ey + exp(1), 1);
        if k>n_discard
            df1dx = 2*x + 2*pi*const_beta*sin(2*pi*x) - const_alpha*ex*(-0.2)*(0.5*x/ts);
            df1dy = -const_alpha*ex*(-0.2)*(0.5*y/ts);
            df2dx = -const_alpha*ey*(-pi*sin(2*pi*x));
            df2dy = 2*y + 2*pi*const_beta*sin(2*pi*y) - const_alpha*ey*(-pi*sin(2*pi*y));
            J = [df1dx, df1dy; df2dx, df2dy];
            [Q, R] = qr(J*Q); s1 = s1 + log(abs(R(1,1))); s2 = s2 + log(abs(R(2,2)));
        end
        x=xn; y=yn;
    end
    lle_4(:, i) = [s1/n_iter_lle; s2/n_iter_lle];
end

f = figure('Visible','off', 'Position', [0 0 1000 600]);
plot(alpha_vals_int, bifur_4, 'k.', 'MarkerSize', 0.1); 
title('Diagram Bifurkacyjny - Odwzorowanie 2D-RA');
xlabel('\alpha'); ylabel('x'); 
xlim([min(alpha_vals_int) max(alpha_vals_int)]);
grid on;
saveas(f, '4_2DRA_bifurkacyjny.png'); close(f);

f = figure('Visible','off', 'Position', [0 0 1000 400]);
plot(alpha_vals_int, lle_4(1,:), 'b-', 'LineWidth', 1.2); hold on; 
plot(alpha_vals_int, lle_4(2,:), 'g--', 'LineWidth', 1.2);
axis tight; ylim('auto');
yl = ylim; m = (yl(2)-yl(1))*0.1; if m==0, m=1; end; ylim([yl(1)-m yl(2)+m]);
legend('\lambda_1', '\lambda_2', 'Location', 'best');
title('Wykładniki Lapunowa - Odwzorowanie 2D-RA'); 
xlabel('\alpha'); ylabel('\lambda'); grid on;
saveas(f, '4_2DRA_LLE.png'); close(f);


%% --- 5. TWO-STAGE LOGISTIC MAP ---
fprintf('[5/7] Generowanie poprawnego Two-Stage Logistic Map...\n');

res_x = 800;          
n_discard = 1000;      
n_keep = 500;         
n_iter_lle = 2000;    

% Parameter r (gamma) in range 3 to 4
r_vals = linspace(3, 4, res_x); 
bifur_5 = zeros(n_keep, res_x); 
lle_5 = zeros(1, res_x);

% If pool is not open, open it
if isempty(gcp('nocreate')), parpool; end

parfor i = 1:res_x
    r = r_vals(i);
    x = 0.1;
    
    % TEMPORARY VARIABLE FOR PARFOR
    % Create vector only for this iteration 'i'
    temp_bifur_col = zeros(n_keep, 1);
    
    % 1. Bifurcation Loop (Transient + Recording)
    for k = 1:n_discard + n_keep
        % Symmetric formula (Double Hump)
        if x < 0.5
            x = 4 * r * x * (0.5 - x);
        else
            % Correct symmetric branch
            x = 4 * r * (1 - x) * (x - 0.5);
        end
        
        % Modulo protection
        if x > 1 || x < 0, x = mod(x, 1); end
        
        % Save to temporary variable
        if k > n_discard
            temp_bifur_col(k - n_discard) = x;
        end
    end
    
    % ASSIGNMENT OF WHOLE COLUMN
    bifur_5(:, i) = temp_bifur_col;
    
    % 2. Calculating LLE
    x = 0.1; 
    sum_L = 0;
    for k = 1:n_iter_lle + n_discard
        % Formula
        if x < 0.5
            xn = 4 * r * x * (0.5 - x);
            deriv = 4 * r * (0.5 - 2*x); % Derivative of left side
        else
            xn = 4 * r * (1 - x) * (x - 0.5);
            deriv = 4 * r * (1.5 - 2*x); % Derivative of right side
        end
        
        if xn > 1 || xn < 0, xn = mod(xn, 1); end
        
        if k > n_discard
            sum_L = sum_L + log(max(abs(deriv), 1e-10));
        end
        x = xn;
    end
    lle_5(i) = sum_L / n_iter_lle;
end

save_plots(5, r_vals, bifur_5, lle_5, '5_SLS', 'a', 'Dwuetapowe odwzorowanie logistyczne (SLS)', [0 1.1]);


%% --- 6. 2D HYPERCHAOTIC MAP ---
fprintf('[6/7] Układ 2D Hyperchaotic Map...\n');
b_bif = 30; b_lle = 40; 
a_vals = linspace(0, 60, res_x);
bifur_6 = zeros(n_keep, res_x); 
lle_6 = zeros(2, res_x);

parfor i = 1:res_x
    a = a_vals(i);
    % Bifurcation (b=30)
    x = 0.2; y = 0.3;
    for k = 1:n_discard
        d=y; if abs(d)<1e-10, d=1e-10; end
        arg1 = a*pi*x+b_bif*y; arg2 = b_bif*pi/d+a*x;
        x=(sin(arg1))^2; y=(cos(arg2))^2;
    end
    for k = 1:n_keep
        d=y; if abs(d)<1e-10, d=1e-10; end
        arg1 = a*pi*x+b_bif*y; arg2 = b_bif*pi/d+a*x;
        xn=(sin(arg1))^2; yn=(cos(arg2))^2;
        bifur_6(k, i) = xn;
        x=xn; y=yn;
    end
    
    % LLE (b=40)
    x=0.2; y=0.3; Q=eye(2); s1=0; s2=0;
    for k=1:n_iter_lle+n_discard
        d=y; if abs(d)<1e-10, d=1e-10; end
        arg1 = a*pi*x+b_lle*y; arg2 = b_lle*pi/d+a*x;
        xn=(sin(arg1))^2; yn=(cos(arg2))^2;
        if k>n_discard
            s2a1=sin(2*arg1); s2a2=sin(2*arg2);
            J = [a*pi*s2a1, b_lle*s2a1; -a*s2a2, (b_lle*pi/(d^2))*s2a2];
            [Q, R] = qr(J*Q); s1=s1+log(abs(R(1,1))); s2=s2+log(abs(R(2,2)));
        end
        x=xn; y=yn;
    end
    lle_6(:, i) = [s1/n_iter_lle; s2/n_iter_lle];
end

f=figure('Visible','off', 'Position', [0 0 1000 600]);
plot(a_vals, bifur_6, 'k.', 'MarkerSize', 0.1); 
title('Diagram Bifurkacyjny - Nowy dwuwymiarowy układ hiperchaotyczny');
xlabel('a'); ylabel('x'); xlim([0 60]); ylim([0 1]);
saveas(f, '6_2DHM_bifurkacyjny.png'); close(f);

f=figure('Visible','off', 'Position', [0 0 1000 400]);
plot(a_vals, lle_6(1,:), 'b'); hold on; plot(a_vals, lle_6(2,:), 'g--');
axis tight; ylim('auto'); legend('\lambda_1', '\lambda_2');
title('Wykładniki Lapunowa - Nowy dwuwymiarowy układ hiperchaotyczny'); 
xlabel('a'); ylabel('\lambda'); grid on;
saveas(f, '6_2DHM_LLE.png'); close(f);


%% --- 7. IMPROVED CHIRIKOV MAP ---
fprintf('[7/7] Improved Chirikov Map...\n');
k_vals = linspace(0, 10, res_x); h_param = 2;
bifur_7 = zeros(n_keep, res_x); 
lle_7 = zeros(1, res_x);

parfor i=1:res_x
    kv = k_vals(i);
    x=0.1; y=0.1;
    for k=1:n_discard
        xn = mod(x + kv*sin(y), 2*pi);
        yn = mod(xn + h_param*y, 2*pi);
        x=xn; y=yn;
    end
    for k=1:n_keep
        xn = mod(x + kv*sin(y), 2*pi);
        yn = mod(xn + h_param*y, 2*pi);
        bifur_7(k, i) = yn; 
        x=xn; y=yn;
    end
    
    x=0.1; y=0.1; Q=eye(2); sum_L=0;
    for k=1:n_iter_lle+n_discard
        xn = mod(x + kv*sin(y), 2*pi);
        yn = mod(xn + h_param*y, 2*pi);
        if k>n_discard
            J = [1, kv*cos(y); 1, kv*cos(y)+h_param];
            [Q, R] = qr(J*Q); sum_L = sum_L + log(abs(R(1,1)));
        end
        x=xn; y=yn;
    end
    lle_7(i) = sum_L/n_iter_lle;
end
save_plots(7, k_vals, bifur_7, lle_7, '7_Chirikov', 'k', 'Ulepszone odwzorowanie Chirikova', [0 2*pi]);

fprintf('\n=== ZAKOŃCZONO. WSZYSTKIE PLIKI ZAPISANE (B&W). ===\n');


%% --- SAVING FUNCTION ---
function save_plots(idx, x_vals, bifur_data, lle, fname, x_lab, title_main, y_lims)
    % Bifurcation Diagram (SCATTER PLOT)
    f = figure('Visible','off', 'Position', [0 0 1000 600]);
    plot(x_vals, bifur_data, 'k.', 'MarkerSize', 0.1);
    title(['Diagram Bifurkacyjny - ' title_main]); 
    xlabel(x_lab); ylabel('x');
    xlim([min(x_vals) max(x_vals)]);
    if ~isempty(y_lims), ylim(y_lims); end
    grid on;
    saveas(f, [fname '_bifurkacyjny.png']); 
    close(f);
    
    % LLE (Standard line plot)
    f = figure('Visible','off', 'Position', [0 0 1000 400]);
    plot(x_vals, lle, 'b-', 'LineWidth', 1.2);
    hold on;
    axis tight; ylim('auto');
    yl = ylim; m = (yl(2)-yl(1))*0.1; if m==0, m=1; end; ylim([yl(1)-m yl(2)+m]);
    title(['Wykładnik Lapunowa - ' title_main]); xlabel(x_lab); ylabel('\lambda'); grid on;
    saveas(f, [fname '_LLE.png']); close(f);
    
    fprintf('Zapisano: %s_*.png\n', fname);
end
