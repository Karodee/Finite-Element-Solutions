%% Free vibration analysis of bar
clc                                                                     %Clears screen
clear all                                                               %Clears workspace

%% Meshing
ne = 2;                                                                 %Number of Elements
nn = 3;                                                                 %Number of nodes
nne = 2;                                                                %Number of nodes per element
dofn = 1;                                                               %Degrees of freedom per node
dofe = 2;                                                               %Degrees of freedom per element
tdof = nn*dofn;                                                         %Degrees of freedom
%CONN = [1,2; 2,3];                                                     %Connectivity Matrix
CONN = zeros(ne, dofe);                                                 %Connectivity Matrix
for ii=1:ne
    %CONN(ii,:) = ii:ii+1;
    CONN(ii,:) = (dofe-1)*(ii-1)+1: dofe + (dofe-1)*(ii-1);
end

%% Geometry and Material properties
le = 2.5;                                                           %Length of elements
Ae = 0.0001;                                                          %Area of elements
E = 210E09;                                                           %Young's Modulus
Rho = 7850;                                                           %Density
P = 5000;                                                             %Concentrated load

%% Assembly of element matrices
KG = zeros(tdof);                                                       %Initializing global stiffness matrix
FGU = zeros(tdof,1);                                                    %Initializing global load matrix  
MG = zeros(tdof);                                                       %Initializing global mass matrix
for i = 1:ne
    KE = E*Ae/le*[1, -1; -1, 1];                                  %Initialize element stiffness matrix
    ME = Rho*Ae*le/6*[2, 1; 1, 2];                                %Initialize element mass matrix
    FE = 0*le/2*[1; 1];
    for j = 1:dofe
        for k = 1:dofe
            KG(CONN(i,j), CONN(i,k)) = KG(CONN(i,j), CONN(i,k)) + KE(j,k);
            MG(CONN(i,j), CONN(i,k)) = MG(CONN(i,j), CONN(i,k)) + ME(j,k);
        end
        FGU(CONN(i,j),1) = FGU(CONN(i,j),1) + FE(j,1);
    end
end
FGC = zeros(tdof,1);
FGC(3,1) = P;
FG = FGU+FGC;
KG = KG+zeros(tdof);

%% Preparing matrices for reaction force analysis
FGR = FG;
KGR = KG;
UGR = zeros(tdof,1);
MGR = MG;

%% Application of Boundary conditions
for kk = 1
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
    MG(kk,:) = [];
    MG(:,kk) = [];
end

%% Solving critical time period
Ktilda = MG\KG;
lambda = eig(Ktilda);
T = 2*pi/(max(lambda))^(1/2);
tcr = T/pi;

%% Initial conditions
delt = 0.00025;
time = 0:0.00025:0.0025;

u=[0;0];
udot=[0;0];
uddot=[0;0];

acc=zeros(11,0);
vel=zeros(11,0);
dis=zeros(11,0);
count=1;


%% Solving Constants
a0 = delt^(-2);
a1 = 1/(2*delt);
a2 = 2*a0;
a3 = 1/a2;

%% Solving
for t = time
    F = 10000-50000*t;
    dis(count) = u(2);
    vel(count) = udot(2);
    acc(count) = uddot(2);
    umindelt = u-delt*udot+a3*uddot;
    Mcap = a0*MG;
    Kcap = KG-a2*MG;
    Rcap = FG - Kcap*u - Mcap*umindelt;
    u = Mcap\Rcap;
    if count>1
        udot = a1*(u-dis(count-1));
        uddot = a0*(dis(count-1)-2*dis(count)+u);
    end
    count=count+1;
end

%% Output
disp('The critical time step for central difference scheme for the given setup in seconds is')
disp(tcr)
saveas(figure,'displacement plot','png')
plot(time,dis)
saveas(figure,'velocity plot','png')
plot(time,vel)
saveas(figure,'accel;eration plot','png')
plot(time,acc)
for i=1:11
    disp('At a time of')
    disp(time(i))
    disp('The displacement, Velocity and Acceleration respectively are as follows')
    disp(dis(i))
    disp(vel(i))
    disp(acc(i))
end
