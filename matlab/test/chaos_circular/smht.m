function out = smht(x,a,b,c,d)
    % SMHT - Helper function for the chaotic map
        out=(exp(a*x)-exp(-b*x))./(exp(c*x)+exp(-d*x));
    end