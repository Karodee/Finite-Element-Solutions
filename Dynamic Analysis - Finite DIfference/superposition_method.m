%% Superposition method
clc
clear all

%% Solving for time
K = [6,-2; -2,4];
M = [2,0; 0,1]; 
Ktilda = M\K;
lambda = eig(Ktilda);
T = 2*pi/(max(lambda))^(1/2);
tcr = T/pi;
delt = T/10;

%% Solving
phi1 = (Ktilda-lambda(1)*[1,0;0,1])\[0;0];
phi2 = (Ktilda-lambda(2)*[1,0;0,1])\[0;0];
OMEGA2 = [lambda(1),0; 0,lambda(2)];
P = zeros(2);
P(1,:)=phi1;
P(2,:)=phi2;