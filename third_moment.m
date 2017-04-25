function [ skew ] = third_moment( data )
%Calculated the skewness for a vector, data

N = length(data);

skew = ((1/(N))*sum((data-mean(data,'omitnan')).^3))/((1/(N-1))*sum((data-mean(data,'omitnan')).^2))^(3/2);


end

