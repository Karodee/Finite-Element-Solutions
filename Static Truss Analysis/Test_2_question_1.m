%% Static analysis of truss using FEM. Assignment 1_1
clc                                     %Clears command window
clear                                   %Clears workspace

%% Matrial properties
E = 200e09;                             %Young's modulus
P1 = -120e03;                           %concentrated load
P2 = -180e03;                            %concentrated load

%% Meshing
ne = 5;                                 %Number of elements
nne = 2;                                %Number of nodes per element
nn = 4;                                 %Number of nodes
dofe = 4;                               %Degrees of freedom per element
dofn = 2;                               %Degrees of freedom per node
tdof = dofn*nn;                         %Total degrees of freedom
CONN = [1 2 3 4; 3 4 5 6; 5 6 7 8; 3 4 7 8; 1 2 7 8];       %Connectivity Matrix
NCONN = zeros(ne, nne);                 %node connectivity matrix
for ii = 1:ne
    NCONN(ii,1) = CONN(ii,2)/2;
    NCONN(ii,2) = CONN(ii,4)/2;
end

%% Geometry
A = 30e-04;                             %uniform area of cross section
coord = [0,0 ; 6,0 ; 10.5,0; 6,6];      %Coordinates of nodes
angles = zeros(ne,1);                   %Angles matrix
le = zeros(ne,1);                       %length of elements
for qq = 1:ne
    angles(qq,1) = atan2((coord(NCONN(qq,2),2)-coord(NCONN(qq,1),2)),(coord(NCONN(qq,2),1)-coord(NCONN(qq,1),1)));
    le(qq,1) = ((coord(NCONN(qq,2),2)-coord(NCONN(qq,1),2))^2 + (coord(NCONN(qq,2),1)-coord(NCONN(qq,1),1))^2)^(1/2);
end

%% Assembly of element matrices
KG = zeros(tdof);                       %Initializing global stiffness matrix
for i = 1:ne
    KE = E*A/le(i)*[1, -1; -1, 1];
    T = [cos(angles(i,1)), sin(angles(i,1)), 0, 0; 0, 0, cos(angles(i,1)), sin(angles(i,1))];
    KEG = T'*KE*T;
    for j = 1:dofe
        for k = 1:dofe
            KG(CONN(i,j), CONN(i,k)) = KG(CONN(i,j), CONN(i,k)) + KEG(j,k);
        end
    end
end
FG = zeros(tdof,1);
FG(7,1) = P1;
FG(8,1) = P2;

%% Preparing matrices for reaction force analysis
FGR = FG;
KGR = KG;
UGR = zeros(tdof,1);

%% Application of boundary conditions
for kk=[1;2;6]
    %disp(kk)
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

%% Final solving
UG = linsolve(KG, FG);
for jj = [3;4;5]
    UGR(jj,1) = UG(jj-2,1);
end
for mm = [7;8]
    UGR(mm,1) = UG(mm-3,1);
end
RF = KGR*UGR-FGR;
BD1 = (UGR(5,1)-UGR(3,1))/4.5 + (UGR(3,1)-UGR(1,1))/6;
BD2 = (UGR(8,1)-UGR(4,1))/6;
BD = E*(BD1^2+BD2^2)^(1/2)*A;

%% Displaying output

disp('Global stifness matrix is')
disp(KGR)
disp('Displacement vector is')
disp(UGR)
disp('Reaction forces at A are')
disp(RF(1:2,1))
disp('Compressive force in BD is')
disp(BD)