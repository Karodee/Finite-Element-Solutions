%% Static analysis of bar using FEM
clc                             %clears screeen
clear all                       %clears workspace


%% Geometry 
ne = 4;                                                                 %Number of elements 
nn = 5;                                                                 %Number of nodes
nne = 3;                                                                %Number of nodes per element
dofn = 2;                                                               %Degrees of freedom per node
dofe = 6;                                                               %Degrees of freedom per element
tdof = dofn*nn;                                                         %Total degrees of freedom
coord = [0,0; 0.4,0; 0.4,0.4; 0,0.4; 0.2,0.2];                          %coordinate matrixAe = 0.04;                                                              %Area of element
h = 0.005;                                                              %Thickness of the plate
Ae = 0.04;
CONN = [1,2,3,4,9,10; 3,4,5,6,9,10; 5,6,7,8,9,10; 7,8,1,2,9,10];        %Connectivity matrix
NCONN = [1,2,5; 2,3,5; 3,4,5; 4,1,5];                                   %Nodal Connectivity Matrix


%% Material Properties and loading
E = 105E09;                                                             %Young's modulus
nu = 0.3;                                                               %Poisson's ratio
D = E/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];                    %D Matrix
Tx = -5E5;                                                               %Traction loading


%% Assembly of element matrices
KG= zeros(tdof, tdof); 
FGU= zeros(tdof,1);
for i=1:ne
    x13 = coord(NCONN(i,1), 1) - coord(NCONN(i,3), 1); 
    x21 = coord(NCONN(i,2), 1) - coord(NCONN(i,1), 1); 
    x32 = coord(NCONN(i,3), 1) - coord(NCONN(i,2), 1); 
    y31 = coord(NCONN(i,3), 2) - coord(NCONN(i,1), 2);
    y12 = coord(NCONN(i,1), 2) - coord(NCONN(i,2), 2);
    y23 = coord(NCONN(i,2), 2) - coord(NCONN(i,3), 2);
    y13 = coord(NCONN(i,1), 2) - coord(NCONN(i,3), 2);
    x23 = coord(NCONN(i,2), 2) - coord(NCONN(i,3), 2); 
    J = 2*(-x21*y23 - x23*y13);                           %Modulus of Jacobian
    Be = [y23, 0, y31, 0, y12, 0; 0, x32, 0, x13, 0, x21; x32, y23, x13, y31, x21, y12]/J;
    Ke = Ae*h*Be'*D*Be;
    for j= 1:dofe
        for k= 1:dofe
            KG(CONN(i,j), CONN(i,k))= KG(CONN(i,j), CONN(i,k))+Ke(j,k);
        end
        %FGU(CONN(i,j),1)= FGU(CONN(i,j),1)+Fe(j,1) + FET(j,1);
    end
end
FGT= zeros(tdof,1);
FGT(3, 1)= Tx/2;
FGT(5, 1)= Tx/2;
FG= FGU+(h*0.4)*FGT;
KGR = KG;
FGR = FG;

%% Application of boundary conditions
for kk=[1;2;7;8]
    %disp(kk)
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

%% Post Pocessing
UG= KG\FG;
UGD = zeros(nn,1); 
UGR = zeros(tdof,1);
UGR(3:6,1) = UG(1:4,1);
UGR(9:10,1) = UG(5:6,1);
stress = zeros(ne,dofe);
for p = 1:nn
    UGD(p,1) = ((UGR(2*p-1,1))^2+(UGR(2*p,1))^2)^(1/2);
end
for r = 1:ne
    Ue = zeros(dofe,1);
    for m = 1:dofe
        Ue(m,1) = UGR(CONN(r,m),1);
    end
    %Ue
    disp('The stress in element')
    disp(r)
    disp(D*Be*Ue)
end

%% Output
disp('The Nodal displacements in x and y direction respectively are')
disp(UGR)
disp('The Nodal displacement magnitudes are')
disp(UGD)
%disp('The stresses in each element are shown in the following Matrix')
%disp(stress)