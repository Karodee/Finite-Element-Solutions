%% Static analysis of bar using FEM
clc                             %clears screeen
clear all                       %clears workspace


%% Geometry 
ne = 12;                                                                 %Number of elements 
nn = 12;                                                                 %Number of nodes
nne = 3;                                                                %Number of nodes per element
dofn = 2;                                                               %Degrees of freedom per node
dofe = 6;                                                               %Degrees of freedom per element
tdof = dofn*nn;                                                         %Total degrees of freedom
%coord = [0,0; 0.4,0; 0.4,0.4; 0,0.4; 0.2,0.2];                          %coordinate matrixAe = 0.04;                                                              %Area of element
h = 0.01;                                                              %Thickness of the plate
Ae = 1.5E-04;
CONN = [1,2,3,4,12,13; 5,6,7,8,15,16; 7,8,9,10,17,18; 1,2,11,12,13,14; 3,4,13,14,15,16; 7,8,15,16,17,18; 11,12,13,14,21,22; 13,14,15,16,23,24; 15,16,17,18,25,26; 11,12,19,20,21,22; 13,14,21,23,23,24; 15,16,23,24,25,26];        %Connectivity matrix
NCONN = [1,2,6; 2,3,7; 3,4,8; 1,5,6; 2,6,7; 3,7,8; 5,6,10; 6,7,11; 7,8,12; 5,9,10; 6,10,11; 7,11,12];                                   %Nodal Connectivity Matrix


%% Material Properties and loading
E = 105E09;                                                             %Young's modulus
nu = 0.3;                                                               %Poisson's ratio
D = E/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];                    %D Matrix
P =-20E3;                                                               %Traction loading


%% Assembly of element matrices
KG= zeros(tdof, tdof); 
FGU= zeros(tdof,1);
for i=1:3
    x13 = 0.02; 
    x21 = 0.02; 
    x32 = 0; 
    y31 = -0.015;
    y12 = 0;
    y23 = -0.015;
    y13 = 0.015;
    x23 = 0; 
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
for i=4:6
    x13 = 0.02; 
    x21 = 0; 
    x32 = 0.02; 
    y31 = -0.015;
    y12 = -0.015;
    y23 = 0;
    y13 = 0.015;
    x23 = 0.02; 
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
for i=7:9
    x13 = 0.02; 
    x21 = 0.02; 
    x32 = 0; 
    y31 = -0.015;
    y12 = 0;
    y23 = -0.015;
    y13 = 0.015;
    x23 = 0; 
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
for i=10:12
    x13 = 0.02; 
    x21 = 0; 
    x32 = 0.02; 
    y31 = -0.015;
    y12 = -0.015;
    y23 = 0;
    y13 = 0.015;
    x23 = 0.02; 
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
FGC= zeros(tdof,1);
FGC(26, 1)= P;
FG= FGU+(h*0.4)*FGC;
KGR = KG;
FGR = FG;

%% Application of boundary conditions
for kk=[1;2;11;12;19;20]
    %disp(kk)
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

UG = KG\FG;
disp('displacement of the point of load')
disp(UG(20,1))