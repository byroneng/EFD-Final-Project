function [ mean ] = first_moment( x )
%calc mean

sum=0;
for i=1:length(x)
    sum=sum+x(i);
end
mean=sum/length(x);
end

