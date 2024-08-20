%% Static analysis of bar using FEM. Assignment 1_1
clc                                     %Clears command window
clear                                   %Clears workspace

%% Matrial properties
E = 200e09;                             %Young's modulus
P1 = -10e03;                           %concentrated load
P2 = 15e03;                            %concentrated load
P3 = -20e03;
P4 = -10e03;

%% Meshing
ne = 10;                                 %Number of elements
nne = 2;                                %Number of nodes per element
nn = 7;                                 %Number of nodes
dofe = 4;                               %Degrees of freedom per element
dofn = 2;                               %Degrees of freedom per node
tdof = dofn*nn;                         %Total degrees of freedom

%% Connectivity
CONN = zeros(ne,dofe);                  %Connectivity Matrix
for ii=1:5
    CONN(ii,:) = 2*ii-1:2*ii+2;
end
CONN(6,1:4) = [11,12,5,6];
CONN(7,1:4) = [5,6,13,14];
CONN(8,1:4) = [13,14,1,2];
CONN(9,1:4) = [3,4,13,14];
CONN(10,1:4) = [7,8,11,12];
NCONN = zeros(ne, nne);                 %node connectivity matrix
for ii = 1:ne
    NCONN(ii,1) = CONN(ii,2)/2;
    NCONN(ii,2) = CONN(ii,4)/2;
end

%% Geometry
A = 3000e-06;                             %uniform area of cross section
coord = [0,-8;4,-4;8,0;12,-3;16,-6;16,0;0,0] ;   %Coordinates of nodes
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
FG(14,1) = P1;
FG(13,1) = P2;
FG(6,1) = P3;
FG(12,1) = P4;

%% Preparing matrices for reaction force analysis
FGR = FG;
KGR = KG;
UGR = zeros(tdof,1);

%% Application of boundary conditions
for kk=[1;2;9;10]
    %disp(kk)
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

%% Final solving
UG = linsolve(KG, FG);
RF = KGR*UGR-FGR;
UGR(3:8,1) = UG(1:6,1);
UGR(11:14,1) = UG(7:10,1);
E1 = (UGR(13,1)-UGR(5,1))/8;
E2 = (UGR(14,1)-UGR(1,1))/8;
EC = -E*(E1^2+E2^2)^(1/2);
G1 = (UGR(11,1)-UGR(5,1))/8;
G2 = (UGR(12,1)-UGR(10,1))/8;
GD = -E*(G1^2+G2^2)^(1/2);

%% Displaying
disp('X and Y displacement at E are')
disp(UG(9:10,1))
disp('X and Y displacement at F are')
disp(UG(3:4,1))
disp('X and Y displacement at G are')
disp(UG(7:8,1))
disp('Compressive stress in EC is')
disp(EC*-1)
disp ('Compressive stress in GD is')
disp(GD*-1)