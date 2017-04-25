function [ kurtosis ] = fourth_moment( x )
%Kurtosis
m=mean(x,'omitnan');
stndDev=second_moment(x)^2;
L=length(x);
% coeff=((L+1)*L)/...
%     ((L-1)*(L-2)*(L-3));
% minus=3*((L-1)^2)/((L-2)*(L-3));
sum=0;

for i=1:L
    sum=(x(i)-m)^4+sum;
end
 kurtosis=(1/L)*(sum/stndDev)-3;


end

