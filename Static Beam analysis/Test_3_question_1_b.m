%% Static analysis of beam using FEM. Assignment 1_1
clc                                     %Clears command window
clear                                   %Clears workspace

%% Matrial properties
E = 200e09;                             %Young's modulus
q0 = zeros(15,1);                       %Distributed load Matrix
for q = 1:10
    q0(q,1) = -30e03;
end
for q = 11:15
    q0(q,1) = 0;
end
P = -55e03;                             %Concentrated load

%% Meshing
ne = 15;                                 %Number of elements
nne = 2;                                %Number of nodes per element
nn = 16;                                 %Number of nodes
dofe = 4;                               %Degrees of freedom per element
dofn = 2;                               %Degrees of freedom per node
tdof = dofn*nn;                         %Total degrees of freedom
CONN = zeros(ne,dofe);                  %Connectivity Matrix
for ii = 1:ne
    CONN(ii,:) = 2*ii-1:2*ii+2;
end
%disp(CONN)

%% Geometry
I = 500e-06;                              %Area moment of inertia
le = 12/15;                                %Length of element

%% Assembly of element matrices
KG = zeros(tdof);                       %Initializing global stiffness matrix
FGU = zeros(tdof,1);                    %Initializing Uniform load matrix
for i = 1:ne
    KE = E*I/(le)^3*[12, 6*le, -12, 6*le; 6*le, 4*(le)^2, -6*le, 2*(le)^2; -12, -6*le, 12, -6*le; -12, 2*(le)^2, -6*(le), 4*(le)^2];
    FE = q0(i)*le/2*[1; le/6; 1; -le/6];
    for j = 1:dofe
        for k = 1:dofe
            KG(CONN(i,j), CONN(i,k)) = KG(CONN(i,j), CONN(i,k)) + KE(j,k);
        end
        FGU(CONN(i,j),1) = FGU(CONN(i,j),1) + FE(j,1);
    end
end
FGC = zeros(tdof,1);                    %Initializing concentrated loads matrix
FGC(31,1) = P;                          %Concentrated load on final node
FG = FGU+FGC;                          %Final load matrix

%% Preparing matrices for reaction force analysis
FGR = FG;
KGR = KG;
UGR = zeros(tdof,1);

%% Application of Boundary conditions
for kk=[1;21]
    %disp(kk)
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

%% Matrix solving
UG = linsolve(KG,FG);                    %Final Displacement Matrix
UGR(2:20,:) = UG(1:19,:);
UGR(22:32,:) = UG(20:30,:);

%% Stress calculation
A = eye(4)\[1, -1, 1, -1; 0, (2/le), (-4/le), (6/le); 1, 1, 1, 1; 0, (2/le), (4/le), (6/le)];
x = zeros(nn,1);
for cc =1:nn
    x(cc) = (cc-1)*le;
end
stress = zeros(nn,1);
for m = 1:15
    We = UGR((2*m-1):(2*m+2), 1);
    strain = 0.1*4/le^2*[0, 0, 2, -6]*A*We;
    stress(m,1) = E*strain;
end
We = UGR(29:32, 1);
strain = 0.1*4/le^2*[0, 0, 2, 6]*A*We;
stress(16,1) = E*strain;
disp(stress)

plot(x, stress)