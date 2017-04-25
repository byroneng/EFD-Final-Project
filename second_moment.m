function [ variance ] = second_moment( data )
%Calculates the variance for the vector data

N = length(data);

variance = (1/(N-1))*sum((data-mean(data,'omitnan')).^2);

end

