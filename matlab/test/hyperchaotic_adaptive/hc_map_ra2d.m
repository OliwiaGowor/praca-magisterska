function seq = hc_map_ra2d(alpha, beta, x, y, len)
    % HC_MAP_RA2D Implementacja mapy chaotycznej 2D Rastrigin-Ackley
    
    bias = 1e8; 
    seq = zeros(len, 1);
    
    % Warm-up
    for i=1:1000
        [x, y] = iteration_ra(x, y, alpha, beta, bias);
    end
    
    % Generowanie właściwe
    count = 0;
    while count < len
        [x, y] = iteration_ra(x, y, alpha, beta, bias);
        
        count = count + 1;
        if count <= len, seq(count) = x; end
        
        count = count + 1;
        if count <= len, seq(count) = y; end
    end
end

function [nx, ny] = iteration_ra(x, y, alpha, beta, bias)
    % Równanie (5) z artykułu
    term1_x = x^2 - (beta + bias) * cos(2*pi*x);
    term2_x = (alpha + bias) * exp(-0.2 * sqrt(0.5*x^2 + 0.5*y^2));
    nx = mod(term1_x - term2_x, 1);
    
    term1_y = y^2 - (beta + bias) * cos(2*pi*y);
    term2_y = (alpha + bias) * exp(0.5*(cos(2*pi*x) + cos(2*pi*y))) + exp(1);
    ny = mod(term1_y - term2_y, 1);
    
    if nx < 0, nx = nx + 1; end
    if ny < 0, ny = ny + 1; end
end