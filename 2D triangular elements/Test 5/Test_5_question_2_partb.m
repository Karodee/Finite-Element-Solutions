%% Static analysis of beam using FEM. Assignment 1_1
clc                                     %Clears command window
clear                                   %Clears workspace

%% Matrial properties
E = 70e09;                             %Young's modulus
%q0 = [0;-12e03];                        %Distributed load

%% Meshing
ne = 10;                                 %Number of elements
nne = 2;                                %Number of nodes per element
nn = 11;                                 %Number of nodes
dofe = 4;                               %Degrees of freedom per element
dofn = 2;                               %Degrees of freedom per node
tdof = dofn*nn;                         %Total degrees of freedom
CONN = zeros(ne,dofe);                  %Connectivity Matrix
for ii = 1:ne
    CONN(ii,:) = 2*ii-1:2*ii+2;
end
%disp(CONN)

%% Geometry
I = 2.25E-08;                              %Area moment of inertia
le = 0.06/10;                                 %Length of element

%% Assembly of element matrices
KG = zeros(tdof);                       %Initializing global stiffness matrix
FGU = zeros(tdof,1);                    %Initializing Uniform load matrix
for i = 1:ne
    KE = E*I/(le)^3*[12, 6*le, -12, 6*le; 6*le, 4*(le)^2, -6*le, 2*(le)^2; -12, -6*le, 12, -6*le; -12, 2*(le)^2, -6*(le), 4*(le)^2];
    %FE = q0(i)*le/2*[1; le/6; 1; -le/6];
    for j = 1:dofe
        for k = 1:dofe
            KG(CONN(i,j), CONN(i,k)) = KG(CONN(i,j), CONN(i,k)) + KE(j,k);
        end
        %FGU(CONN(i,j),1) = FGU(CONN(i,j),1) + FE(j,1);
    end
end
FGC = zeros(tdof,1);                    %Initializing concentrated loads matrix
FGC(22,1) = -2E5;
FG = FGU+FGC;                           %Final load matrix

%% Preparing matrices for reaction force analysis
FGR = FG;
KGR = KG;
UGR = zeros(tdof,1);

%% Application of Boundary conditions
for kk=[1;2]
    %disp(kk)
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

%% Solving
UG= KG\FG;

%%Post Procsessing
disp('The displacement of node of contact with the concentrated load is')
disp(((UG(19,1))^2+(UG(20,1))^2)^(1/2))