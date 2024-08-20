%% Static analysis of frames using FEM. Assignment 1_1
clc                                     %Clears command window
clear                                   %Clears workspace

%% Matrial properties
E = 200e09;                             %Young's modulus
q0 = [0;0;-20e03];                        %Distributed load
P = -150e03;

%% Geometry
I = 500e-06;                              %Area moment of inertia
coord = [0,0 ; 1.5,2 ; 3,4; 8,4];         %Coodsinates of nodes
A = 0.15;                                 %Assuming the beam elements are rectangular

%% Meshing
ne = 3;                                 %Number of elements
nne = 2;                                %Number of nodes per element
nn = 4;                                 %Number of nodes
dofe = 6;                               %Degrees of freedom per element
dofn = 3;                               %Degrees of freedom per node
tdof = dofn*nn;                         %Total degrees of freedom
CONN = zeros(ne,dofe);                  %Connectivity Matrix
for ii = 1:ne
    CONN(ii,:) = 2*ii-1:2*ii+4;
end
NCONN = zeros(ne, nne);                 %node connectivity matrix
for ii = 1:ne
    NCONN(ii,1) = CONN(ii,2)/2;
    NCONN(ii,2) = CONN(ii,4)/2;
end
angles = zeros(ne,1);                   %Angles matrix
le = zeros(ne,1);                       %length of elements
for qq = 1:ne
    angles(qq,1) = atan2((coord(NCONN(qq,2),2)-coord(NCONN(qq,1),2)),(coord(NCONN(qq,2),1)-coord(NCONN(qq,1),1)));
    le(qq,1) = ((coord(NCONN(qq,2),2)-coord(NCONN(qq,1),2))^2 + (coord(NCONN(qq,2),1)-coord(NCONN(qq,1),1))^2)^(1/2);
end
%disp(CONN)

%% Assembly of element matrices
KG = zeros(tdof);                       %Initializing global stiffness matrix
FGU = zeros(tdof,1);                    %Initializing Uniform load matrix
for i = 1:ne
    KE = [E*A/le(i), 0, 0, -E*A/le(i), 0, 0; 0, 12*E*I/le(i)^3, 6*E*I/le(i)^2, 0, -12*E*I/le(i)^3, 6*E*I/le(i)^2; 0, 6*E*I/le(i)^2, 4*E*I/le(i)^2, 0, -6*E*I/le(i)^2, 2*E*I/le(i); -E*A/le(i), 0, 0, E*A/le(i), 0, 0; 0, -12*E*I/le(i)^3, -6*E*I/le(i)^2, 0, 12*E*I/le(i)^3, -6*E*I/le(i)^2; 0, 6*E*I/le(i)^2, 2*E*I/le(i)^2, 0, -6*E*I/le(i)^2, 4*E*I/le(i)];
    FE = q0(i)*le(i)/2*[0; 1; le(i)/6;0;  1; -le(i)/6];
    T = [cos(angles(i,1)), sin(angles(i,1)), 0, 0, 0, 0; -sin(angles(i,1)), cos(angles(i,1)), 0, 0, 0, 0; 0, 0, 1, 0, 0, 0; 0, 0, 0, cos(angles(i,1)), sin(angles(i,1)), 0; 0, 0, 0, -sin(angles(i,1)), cos(angles(i,1)), 0; 0, 0, 0, 0, 0, 1];
    KEG = T'*KE*T;
    for j = 1:dofe
        for k = 1:dofe
            KG(CONN(i,j), CONN(i,k)) = KG(CONN(i,j), CONN(i,k)) + KE(j,k);
        end
        FGU(CONN(i,j),1) = FGU(CONN(i,j),1) + FE(j,1);
    end
end
FGC = zeros(tdof,1);                    %Initializing concentrated loads matrix
FGC(3,1) = P;
FG = FGU+FGC;                           %Final load matrix

%% Preparing matrices for reaction force analysis
FGR = FG;
KGR = KG;
UGR = zeros(tdof,1);

%% Application of Boundary conditions
for kk=[1;2;11]
    %disp(kk)
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

%% Matrix solving
UG = linsolve(KG,FG)                    %Final Displacement Matrix
UGR(4,:) = UG(1,:);
UGR(6,:) = UG(2,:);
RF = KGR*UGR-FGR                        %Reaction Forces Matrix
