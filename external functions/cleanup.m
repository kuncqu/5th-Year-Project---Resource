function [ N, F, A ] = cleanup( N, F, A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


for i = 1:length(N)
     if norm(N(i,:)) == 0 || A(i,1) == 0  
        N(i,:) = [0 0 0];
        F(i,:) = [0 0 0];
        A(i,:) = 0;    
     end     
end

N( ~any(N,2), : ) = [];
F( ~any(F,2), : ) = [];
A( ~any(A,2), : ) = [];

end

