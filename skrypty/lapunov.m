clear; clc; close all;

%% === USTAWIENIA GLOBALNE ===
res_x = 600;          % Rozdzielczość pozioma (dla 2D-RA to liczba całkowitych kroków)
res_y = 800;          % Rozdzielczość pionowa (jakość diagramu gęstości)
n_iter_bif = 5000;    % Liczba iteracji dla diagramu (im więcej tym lepszy kontrast)
n_iter_lle = 1000;    % Liczba iteracji dla LLE
n_discard = 500;      % Odrzucenie stanu przejściowego

% Konfiguracja puli wątków dla parfor (opcjonalne, przyspiesza obliczenia)
if isempty(gcp('nocreate')), parpool; end

fprintf('=== START SYMULACJI (7 UKŁADÓW) ===\n');

%% --- 1. UKŁAD KOSINUSOWO-WIELOMIANOWY ---
fprintf('[1/7] Układ Kosinusowo-Wielomianowy...\n');
a_cos = 4;
mu_vals = linspace(0, 5, res_x);
y_lims = [-1.1, 1.1]; 
y_edges = linspace(y_lims(1), y_lims(2), res_y + 1);

dens_1 = zeros(res_y, res_x);
lle_1 = zeros(1, res_x);

parfor i = 1:res_x
    mu = mu_vals(i);
    % Bifurkacja
    x = 0.1; orbit = zeros(n_iter_bif, 1);
    for k = 1:n_iter_bif+n_discard
        arg = mu*(x^3+x) + a_cos; x = cos(arg);
        if k>n_discard, orbit(k-n_discard)=x; end
    end
    dens_1(:, i) = histcounts(orbit, y_edges)';
    % LLE
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
save_plots(1, mu_vals, y_lims, dens_1, lle_1, '1_CosPoly', '\mu', '1D Cos-Poly', false);


%% --- 2. 3D LOGISTIC MAP ---
fprintf('[2/7] Układ 3D Logistic Map...\n');
beta = 0.017; gamma = 0.0012;
alpha_vals = linspace(3.65, 4.0, res_x);
y_lims = [0, 1];
y_edges = linspace(0, 1, res_y+1);

dens_2 = zeros(res_y, res_x);
lle_2 = zeros(1, res_x);

parfor i = 1:res_x
    alpha = alpha_vals(i);
    x=0.2487; y=0.35; z=0.10001; orbit=zeros(n_iter_bif,1);
    % Bifurkacja
    for k=1:n_iter_bif+n_discard
        xn = alpha*x*(1-x) + beta*y^2*x + gamma*z^3;
        yn = alpha*y*(1-y) + beta*z^2*y + gamma*x^3;
        zn = alpha*z*(1-z) + beta*x^2*z + gamma*y^3;
        if k>n_discard, orbit(k-n_discard)=xn; end
        x=xn; y=yn; z=zn;
    end
    dens_2(:, i) = histcounts(orbit, y_edges)';
    % LLE
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
save_plots(2, alpha_vals, y_lims, dens_2, lle_2, '2_3DLogistic', '\alpha', '3D Logistic', false);


%% --- 3. SOBOLEVA-MODULO MAP ---
fprintf('[3/7] Układ Soboleva-Modulo...\n');
a_sob=5; b_sob=5; K_sob=1; B_sob=1; C_sob=1; D_sob=1;
A_vals = linspace(0, 5, res_x);
y_lims = [0, 1]; y_edges = linspace(0, 1, res_y+1);

dens_3 = zeros(res_y, res_x);
lle_3 = zeros(1, res_x);

parfor i = 1:res_x
    A_val = A_vals(i);
    x=0.1; orbit=zeros(n_iter_bif,1);
    for k=1:n_iter_bif+n_discard
        u = x/K_sob;
        smht = (exp(A_val*u)-exp(-B_sob*u))/(exp(C_sob*u)+exp(-D_sob*u));
        xn = mod(x + a_sob + b_sob*smht, K_sob);
        if k>n_discard, orbit(k-n_discard)=xn; end
        x=xn;
    end
    dens_3(:, i) = histcounts(orbit, y_edges)';
    % LLE
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
save_plots(3, A_vals, y_lims, dens_3, lle_3, '3_Soboleva', 'A', 'Soboleva-Modulo', false);


%% --- 4. 2D-RA MAP (SPECIAL: Integers, Bias, 2 LEs) ---
fprintf('[4/7] Układ 2D-RA Map (Bias 10^8, Parametry Całkowite)...\n');
bias = 1e8; beta_ra = 1;
alpha_vals = 0:(res_x-1); % Liczby całkowite
y_lims = [0, 1]; y_edges = linspace(0, 1, res_y+1);

dens_4 = zeros(res_y, res_x);
lle_4 = zeros(2, res_x); % Dwa wykładniki

parfor i = 1:res_x
    alpha = alpha_vals(i);
    const_beta = beta_ra + bias; const_alpha = alpha + bias;
    x=0.5; y=0.5; orbit=zeros(n_iter_bif,1);
    
    % Bifurkacja
    for k=1:n_iter_bif+n_discard
        ts = sqrt(0.5*x^2+0.5*y^2); if ts<1e-10, ts=1e-10; end
        ex = exp(-0.2*ts); ey = exp(0.5*cos(2*pi*x)+0.5*cos(2*pi*y));
        xn = mod(x^2 - const_beta*cos(2*pi*x) - const_alpha*ex, 1);
        yn = mod(y^2 - const_beta*cos(2*pi*y) - const_alpha*ey + exp(1), 1);
        if k>n_discard, orbit(k-n_discard)=xn; end
        x=xn; y=yn;
    end
    dens_4(:, i) = histcounts(orbit, y_edges)';
    
    % 2x LLE
    x=0.5; y=0.5; Q=eye(2); s1=0; s2=0;
    for k=1:n_iter_lle+n_discard
        ts = sqrt(0.5*x^2+0.5*y^2); if ts<1e-10, ts=1e-10; end
        ex = exp(-0.2*ts); ey = exp(0.5*cos(2*pi*x)+0.5*cos(2*pi*y));
        xn = mod(x^2 - const_beta*cos(2*pi*x) - const_alpha*ex, 1);
        yn = mod(y^2 - const_beta*cos(2*pi*y) - const_alpha*ey + exp(1), 1);
        
        if k>n_discard
            df1dx = 2*x + 2*pi*const_beta*sin(2*pi*x) - const_alpha*ex*(-0.2)*(0.5*x/ts);
            df1dy = -const_alpha*ex*(-0.2)*(0.5*y/ts);
            df2dx = -const_alpha*ey*(-pi*sin(2*pi*x)); % Poprawione pochodne wewn.
            df2dy = 2*y + 2*pi*const_beta*sin(2*pi*y) - const_alpha*ey*(-pi*sin(2*pi*y));
            J = [df1dx, df1dy; df2dx, df2dy];
            [Q, R] = qr(J*Q);
            s1 = s1 + log(abs(R(1,1))); s2 = s2 + log(abs(R(2,2)));
        end
        x=xn; y=yn;
    end
    lle_4(:, i) = [s1/n_iter_lle; s2/n_iter_lle];
end

% Manualny zapis dla 2D-RA (specyficzny wykres)
f=figure('Visible','off'); imagesc(alpha_vals, y_lims, log10(dens_4+1)); 
set(gca,'YDir','normal'); colormap(jet(256)); shading flat; title('2D-RA: Diagram Bifurkacyjny');
xlabel('Parametr \alpha (int)'); ylabel('x_n'); saveas(f, '4_2DRA_bifurkacyjny.png'); close(f);

f=figure('Visible','off'); plot(alpha_vals, lle_4(1,:), 'b'); hold on; plot(alpha_vals, lle_4(2,:), 'g--');
yline(0,'r'); axis tight; ylim('auto'); legend('\lambda_1', '\lambda_2'); title('2D-RA: 2 Wykładniki (Hiperchaos)');
xlabel('Parametr \alpha (int)'); ylabel('\lambda'); grid on; saveas(f, '4_2DRA_LLE.png'); close(f);
fprintf('Zapisano: 4_2DRA_bifurkacyjny.png, 4_2DRA_LLE.png\n');


%% --- 5. 1D SLS MAP ---
fprintf('[5/7] Układ 1D SLS Map...\n');
a_vals = linspace(0, 1, res_x);
y_lims = [-1, 1]; y_edges = linspace(-1, 1, res_y+1);

dens_5 = zeros(res_y, res_x);
lle_5 = zeros(1, res_x);

parfor i=1:res_x
    a = a_vals(i);
    x=0.5; orbit=zeros(n_iter_bif,1);
    for k=1:n_iter_bif+n_discard
        arg = pi*(4*a*x*(1-x) + (1-a)*sin(pi*x));
        xn = sin(arg);
        if k>n_discard, orbit(k-n_discard)=xn; end
        x=xn;
    end
    dens_5(:, i) = histcounts(orbit, y_edges)';
    % LLE
    x=0.5; sum_L=0;
    for k=1:n_iter_lle+n_discard
        arg = pi*(4*a*x*(1-x) + (1-a)*sin(pi*x));
        xn = sin(arg);
        if k>n_discard
            u_prime = pi*(4*a*(1-2*x) + (1-a)*pi*cos(pi*x));
            deriv = cos(arg)*u_prime;
            sum_L = sum_L + log(max(abs(deriv), 1e-10));
        end
        x=xn;
    end
    lle_5(i) = sum_L/n_iter_lle;
end
save_plots(5, a_vals, y_lims, dens_5, lle_5, '5_SLS', 'a', '1D SLS', false);


%% --- 6. 2D HYPERCHAOTIC MAP (Dopasowany do Fig. 3c) ---
fprintf('[6/7] Układ 2D Hyperchaotic Map (b=40, 2 wykładniki)...\n');

% Parametry zgodne z tekstem: "Figure 3c ... parameter b = 40 while a varies between 0 and 60"
b_hm = 40; 

a_vals = linspace(0, 60, res_x);
y_lims = [0, 1]; 
y_edges = linspace(0, 1, res_y+1);

dens_6 = zeros(res_y, res_x);
lle_6 = zeros(2, res_x); % Miejsce na dwa wykładniki

parfor i = 1:res_x
    a = a_vals(i);
    
    % --- 1. Bifurkacja (Gęstość) ---
    x = 0.2; y = 0.3; 
    orbit = zeros(n_iter_bif, 1);
    
    for k = 1:(n_iter_bif + n_discard)
        % Zabezpieczenie mianownika
        denom = y; if abs(denom) < 1e-10, denom = 1e-10; end
        
        arg1 = a * pi * x + b_hm * y;
        arg2 = b_hm * pi / denom + a * x;
        
        % Równania (Eq. 1 z artykułu)
        xn = (sin(arg1))^2;
        yn = (cos(arg2))^2;
        
        if k > n_discard
            orbit(k - n_discard) = xn; 
        end
        x = xn; y = yn;
    end
    dens_6(:, i) = histcounts(orbit, y_edges)';
    
    % --- 2. Dwa Wykładniki Lapunowa ---
    x = 0.2; y = 0.3; 
    Q = eye(2); 
    s1 = 0; s2 = 0;
    
    for k = 1:(n_iter_lle + n_discard)
        denom = y; if abs(denom) < 1e-10, denom = 1e-10; end
        
        arg1 = a * pi * x + b_hm * y;
        arg2 = b_hm * pi / denom + a * x;
        
        xn = (sin(arg1))^2;
        yn = (cos(arg2))^2;
        
        if k > n_discard
            % Obliczanie Macierzy Jacobiego
            sin2_arg1 = sin(2 * arg1);
            sin2_arg2 = sin(2 * arg2);
            
            j11 = a * pi * sin2_arg1;
            j12 = b_hm * sin2_arg1;
            j21 = -a * sin2_arg2;
            j22 = (b_hm * pi / (denom^2)) * sin2_arg2;
            
            J = [j11, j12; j21, j22];
            
            % Dekompozycja QR
            [Q, R] = qr(J * Q);
            s1 = s1 + log(abs(R(1,1)));
            s2 = s2 + log(abs(R(2,2)));
        end
        x = xn; y = yn;
    end
    % Zapisujemy oba wyniki jako wektor kolumnowy
    lle_6(:, i) = [s1 / n_iter_lle; s2 / n_iter_lle];
end

% --- ZAPIS WYKRESÓW DLA UKŁADU 6 ---

% Diagram Bifurkacyjny
f = figure('Visible','off', 'Position', [0 0 1000 600]);
imagesc(a_vals, y_lims, log10(dens_6 + 1));
set(gca, 'YDir', 'normal');
colormap(jet(256)); shading interp; colorbar;
title('2D Hyperchaotic Map: Diagram Bifurkacyjny (b=40)');
xlabel('Parametr a'); ylabel('x_n');
saveas(f, '6_2DHM_bifurkacyjny.png'); 
close(f);

% Wykres Dwóch Wykładników (LLE)
f = figure('Visible','off', 'Position', [0 0 1000 400]);
plot(a_vals, lle_6(1,:), 'b-', 'LineWidth', 1.5, 'DisplayName', '\lambda_1');
hold on;
plot(a_vals, lle_6(2,:), 'g--', 'LineWidth', 1.5, 'DisplayName', '\lambda_2');
yline(0, 'r-', 'LineWidth', 1.0);
legend('show');
axis tight; ylim('auto');
title('2D Hyperchaotic Map: Spektrum Wykładników (b=40)');
xlabel('Parametr a'); ylabel('\lambda'); grid on;
saveas(f, '6_2DHM_LLE.png'); 
close(f);

fprintf('Zapisano: 6_2DHM_bifurkacyjny.png, 6_2DHM_LLE.png\n');

%% --- 7. IMPROVED CHIRIKOV MAP ---
fprintf('[7/7] Improved Chirikov Map...\n');
k_vals = linspace(0, 10, res_x); h_param = 2;
y_lims = [0, 2*pi]; y_edges = linspace(0, 2*pi, res_y+1);
dens_7 = zeros(res_y, res_x); lle_7 = zeros(1, res_x);

parfor i=1:res_x
    kv = k_vals(i);
    x=0.1; y=0.1; orbit=zeros(n_iter_bif,1);
    for k=1:n_iter_bif+n_discard
        xn = mod(x + kv*sin(y), 2*pi);
        yn = mod(xn + h_param*y, 2*pi);
        if k>n_discard, orbit(k-n_discard)=yn; end
        x=xn; y=yn;
    end
    dens_7(:, i) = histcounts(orbit, y_edges)';
    % LLE
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
save_plots(7, k_vals, y_lims, dens_7, lle_7, '7_Chirikov', 'k', 'Improved Chirikov', false);

fprintf('\n=== ZAKOŃCZONO. WSZYSTKIE PLIKI ZAPISANE. ===\n');


%% --- FUNKCJA ZAPISUJĄCA ---
function save_plots(idx, x_vals, y_lims, dens, lle, fname, x_lab, title_main, is_int)
    % Diagram Bifurkacyjny
    f = figure('Visible','off', 'Position', [0 0 1000 600]);
    imagesc(x_vals, y_lims, log10(dens + 1));
    set(gca, 'YDir', 'normal');
    colormap(jet(256)); colorbar; 
    if is_int, shading flat; else, shading interp; end
    title([title_main ': Diagram Bifurkacyjny']); xlabel(['Parametr ' x_lab]); ylabel('Stan');
    saveas(f, [fname '_bifurkacyjny.png']); close(f);
    
    % LLE
    f = figure('Visible','off', 'Position', [0 0 1000 400]);
    plot(x_vals, lle, 'b-', 'LineWidth', 1.2);
    hold on; yline(0, 'r--', 'LineWidth', 1.5);
    axis tight; ylim('auto');
    yl = ylim; m = (yl(2)-yl(1))*0.1; if m==0, m=1; end; ylim([yl(1)-m yl(2)+m]);
    title([title_main ': LLE']); xlabel(['Parametr ' x_lab]); ylabel('\lambda'); grid on;
    saveas(f, [fname '_LLE.png']); close(f);
    
    fprintf('Zapisano: %s_*.png\n', fname);
end