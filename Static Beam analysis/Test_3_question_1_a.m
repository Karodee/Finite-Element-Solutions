%% Static analysis of beam using FEM
clc                                     %Clears command window
clear                                   %Clears workspace

%% Matrial properties
E = 200e09;                             %Young's modulus
q0 = [-30e03; 0];                       %Distributed load
P = -55e03;                             %Concentrated load

%% Meshing
ne = 2;                                 %Number of elements
nne = 2;                                %Number of nodes per element
nn = 3;                                 %Number of nodes
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
le = [9,3];                                %Length of element

%% Assembly of element matrices
KG = zeros(tdof);                       %Initializing global stiffness matrix
FGU = zeros(tdof,1);                    %Initializing Uniform load matrix
for i = 1:ne
    KE = E*I/(le(i))^3*[12, 6*le(i), -12, 6*le(i); 6*le(i), 4*(le(i))^2, -6*le(i), 2*(le(i))^2; -12, -6*le(i), 12, -6*le(i); -12, 2*(le(i))^2, -6*(le(i)), 4*(le(i))^2];
    FE = q0(i)*le(i)/2*[1; le(i)/6; 1; -le(i)/6];
    for j = 1:dofe
        for k = 1:dofe
            KG(CONN(i,j), CONN(i,k)) = KG(CONN(i,j), CONN(i,k)) + KE(j,k);
        end
        FGU(CONN(i,j),1) = FGU(CONN(i,j),1) + FE(j,1);
    end
end
FGC = zeros(tdof,1);                    %Initializing concentrated loads matrix
FGC(5,1) = P;                          %Concentrated load on final node
FG = FGU+FGC;                          %Final load matrix

%% Application of Boundary conditions
for kk=[1;3]
    %disp(kk)
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

%% Matrix solving
UG = linsolve(KG,FG)                    %Final Displacement Matrix

%% Output
disp('Displacement at point C in m is')
disp(UG(3,1))
disp('Slope at point B is')
disp(UG(2,1))