%% Newmark Method
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

%% Solving Constants
alpha=0.25;
delta=0.5;
a0 = 1/(alpha*delt^2);
a1 = delta/(alpha*delt);
a2 = 1/(alpha*delt);
a3 = 1/(2*alpha)-1;
a4 = delta/alpha-1;
a5 = delt/2*(delta/alpha-2);
a6 = delt*(1-delta);
a7 = delta*delt;
