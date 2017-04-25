function[rho] = AutoCorrelation(x)

x=inpaint_nans(x);
variance =0;
top=0;

N=length(x);
x = x - mean(x,'omitnan');

for i =1:N
     variance = variance + (x(i))^2;
end

variance = variance/N;

for k = 0:N/2 %k is the time lag

    for i= 1:(N-k)
        top = top +(x(i))*(x(i+k));
    end
    
    top = top/(length(x)-k);
    rho(k+1) = top/variance;
end
end




