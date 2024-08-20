%% Central Difference Method
clc
clear all

%% Initial conditions
K = 20000;
M = 6000; 
delt = 0.05;
time=[0,0.05,0.1,0.15];
F = 10000-50000*t;

u=0;
udot=0;
uddot=0;

acc=zeros(4,0);
vel=zeros(4,0);
dis=zeros(4,0);
count=0;


%% Solving Constants
a0 = delt^(-2);
a1 = 1/(2*delt);
a2 = 2*a0;
a3 = 1/a2;

%% Solving
for t = time
    dis(count) = u;
    vel(count) = udot;
    acc(count) = uddot;
    umindelt = u-delt*udot+a3*uddot;
    Mcap = a0*M;
    Kcap = K-a2*M;
    Rcap = F - Kcap*u - Mcap*umindelt;
    u = Mcap\Rcap;
    if count>0
        udot = a1*(u-dis(count-1));
        uddot = a0*(dis(count-1)-2*dis(count)+u);
    end
    count=count+1;
end

plot()
