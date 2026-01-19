function seq = two_stage_logistic(x0, gamma, len)
    % Generates a chaotic sequence of length 'len'
    % Formula from Section 3.1.2 
    
    seq = zeros(len, 1);
    x = x0;
    
    for k = 1:len
        if x >= 0 && x < 0.5
            x = 4 * gamma * x * (0.5 - x);
        elseif x >= 0.5 && x <= 1
            x = 1 - 4 * gamma * x * (x - 0.5);
        end
        seq(k) = x;
    end
end
