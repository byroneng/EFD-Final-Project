function [ var ] = rearrangeHeights( var )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

temp=var;
var(:,3)=temp(:,6);
var(:,4)=temp(:,3);
var(:,5)=temp(:,4);
var(:,6)=temp(:,5);
end

