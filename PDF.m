function [ pdf ] = PDF( std, mean, x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

pdf=(std.*(2.*pi).^(.5)).^(-1).*exp((-(x-mean).^2)./(2.*std.^2));


end

