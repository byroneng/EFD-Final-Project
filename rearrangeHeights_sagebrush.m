function [ var ] = rearrangeHeights_sagebrush( var )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

temp=var;
var(:,2)=temp(:,5);
var(:,3)=temp(:,2);
var(:,4)=temp(:,3);
var(:,5)=temp(:,4);
end
