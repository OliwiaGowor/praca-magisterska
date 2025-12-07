function out = smht(x,a,b,c,d)

out=(exp(a*x)-exp(-b*x))./(exp(c*x)+exp(-d*x));
end