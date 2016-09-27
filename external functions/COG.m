function [ CG ] = COG( P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

X = sum(P(:,1))/length(P);
Y = sum(P(:,2))/length(P);
Z = sum(P(:,3))/length(P);

CG = [X,Y,Z];

end

