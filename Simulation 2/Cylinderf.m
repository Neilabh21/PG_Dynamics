function [ dy ] = Cylinderf( t,y )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% These parameters should be actually in a different file and defined as
% global variables to be shared (since this is small simulation we can
% afford to spend time in trivially defining them every time) 
m=0.2;% 0.10; % Kg
l=0.15;  % m
g=9.81;
I0 = 1/12*m*l^2;
F = 0.1; % N

dy = zeros(3,1);
dy(2) = y(3);

dy(3) = F * l / (2 * I0);
dy(1) = F/m;
 end

