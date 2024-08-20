%% Static Analysis of plate using CPT quadrilateral elements
clc                             %clears screeen
clear all                       %clears workspace

%% Meshing
ne = 4;                                         %Number of elements
nne = 4;                                        %Numeber of nodes per element
nn = 9;                                         %Number of nodes
dofn = 5;                                       %Degrees of freedom per node
dofe = 20;                                      %Degrees of freedom per element
tdof = dofn*nn;                                 %Total Degrees of Freedom
NCONN = [1,2,5,4; 2,3,6,5; 4,5,8,7; 5,6,9,8];   %Nodal Connectivity matrix
CONN=zeros(ne,dofe);                            %Connectivity Matrix
for i=1:ne
    for j=1:nne
        for k=1:dofn
            CONN(i,(j-1)*5+k) = (NCONN(i,j)-1)*5+k
        end
    end
end
coord = [0,0; 1,0; 2,0; 0,1; 1,1; 2,1; 0,2; 1,2; 2,2];

%% Geometry and Material Properties
Ae = 0.01;                                      %Area of Element
E = 200E09;                                     %Young's modulus
nu = 0.25;                                      %Poisson's ratio
t = 0.05;                                       %thickness of plate
lx = 1;                                         %length of element in x-direction
ly = 1;                                         %ength of element in y-direction
q0 = [0,0,100];                                 %uniformly distributed load

%% Assembly of element stiffness matrices
% 4 point formula
Gp=[-0.8611 -0.3399 0.3399 0.8611];             %Gauss Points
W=[0.3478 0.652 0.652 0.3478];

p=size(Gp,2);                                   %No. of columns in Gp
q=size(Gp,2);
D=(E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2];
KG=zeros(tdof,tdof);
FG = zeros(tdof,1);
for i=1:ne                                      %Membrane stiffness
    de=[coord(NCONN(i,1),1);coord(NCONN(i,1),2);coord(NCONN(i,2),1);coord(NCONN(i,2),2);coord(NCONN(i,3),1);coord(NCONN(i,3),2);coord(NCONN(i,4),1);coord(NCONN(i,4),2)];
    kem=zeros(dofe,dofe);
    keb=zeros(dofe,dofe);
    Ke=zeros(dofe);
    Fe = zeros(dofe,1);
    for k = 1:p
        for l = 1:q
            e=Gp(k);    % Taking Gauss points into e and n
            n=Gp(l);

            phi = 0.25*[(1-e)(1-n), (1+e)(1-n), (1+e)(1+n), (1-e)(1+n)];
            phi_z = 0.25*[-(1-n), (1-n), (1+n), -(1+n)];
            phi_n = 0.25*[-(1-e), -(1+e), (1+e), (1-e)];
            X = [de(1,1); de(2,1); de(3,1); de(4,1)];
            Y = [de(1,2); de(2,2); de(3,2); de(4,2)];


            %Membrane stiffness matrix
            x_z=0.25*[-(1-n), (1-n), (1+n), -(1+n)]*X;   
            x_n=0.25*[-(1-e), -(1+e), (1+e), (1-e)]*X;
            y_z=0.25*[-(1-n), (1-n), (1+n), -(1+n)]*Y;
            y_n=0.25*[-(1-e), -(1+e), (1+e), (1-e)]*Y;
            
            J=x_z*y_n-x_n*y_z;
            phi_x = y_n/J*phi_z-y_z/J*phi_n;
            phi_y = -x_n/J*phi_z+x_z/J*phi_n;
            Bme = zeros(3,dofe);
            for a=1:nne
                Bme(1,a) = phi_x(1,a);
                Bme(2,a+4) = phi_y(1,a);
                Bme(3,a) = phi_x(1,a);
                Bme(3,a+4) = phi_y(1,a);
            end
            kem =t*Bme'*D*Bme*J*W(k)*W(l)+kem;                      %Membrane element stiffness matrix

            %bending stiffness matrix
            A = [1,-1,-1,1,1,1,-1,-1,-1,-1,1,1; 1,1,-1,1,-1,1,1,-1,1,-1,-1,-1; 1,1,1,1,1,1,1,1,1,1,1,1; 1,-1,1,1,-1,1,-1,1,-1,1,-1,-1; 0,1/lx,0,-2/lx,-1/lx,0,3/lx,2/lx,1/lx,0,-3/lx,-1/lx; 0,1/lx,0,2/lx,-1/lx,0,3/lx,-2/lx,1/lx,0,-3/lx,-1/lx; 0,1/lx,0,2/lx,1/lx,0,3/lx,2/lx,1/lx,0,3/lx,1/lx; 0,1/lx,0,-2/lx,1/lx,0,3/lx,-2/lx,1/lx,0,3/lx,1/lx; 0,0,1/ly,0,-1/ly,-2/ly,0,1/ly,2/ly,3/ly,-1/ly,-3/ly; 0,0,1/ly,0,1/ly,-2/ly,0,1/ly,-2/ly,3/ly,1/ly,3/ly; 0,0,1/ly,0, 1/ly,2/ly,0,1/ly,-2/ly,3/ly,1/ly,3/ly; 0,0,1/ly,0,-1/ly,2/ly,0,1/ly,-2/ly,3/ly,-1/ly,-3/ly];
            Psi = [1,e,n,e^2,e*n,n^2,e^3,e^2*n,e*n^2,n^3,e^3*n,e*n^3];
            Psi_z = [0,1,0,2e,n,0,3*e^2,2*e*n,n^2,0,3*e^2*n,n^3];
            Psi_n = [0,0,1,0,e,2*n,0,2*e^2,2*e*n,3*n^2,e^3,3*n^2*e];
            Psi_zz = [0,0,0,2,0,0,6*e,2*n,0,0,6*e*n,0];
            Psi_nn = [0,0,0,0,0,2,0,0,2*e,6*n,0,6*n*e];
            
            B = Psi*inv(A);                                         %defining Beta matrices
            Jw = [x_z,y_z,0,0,0; x_n,y_n,0,0,0; 0,0,x_z^2,y_z^2,2*x_z*y_z; 0,0,x_n^2,x_y^2,2*x_n*y_n; 0,0,x_n*x_z,y_n*y_z,x_z*y_n+x_n*y_z];
            Psicart = inv(Jw)*[Psi;Psi_z;Psi_n;Psi_zz;Psi_nn];
            Bxx = Psicart(3,1:12)*inv(A);
            Byy = Psicart(4,1:12)*inv(A);
            Bxy = Psicart(5,1:12)*inv(A);
            
            Bbe = zeros(3,20);                                  
            Bbe(1,9:20)=Bxx;
            Bbe(2,9:20)=Byy;
            Bbe(3,9:20)=Bxy;
            
            keb = t^3/12*Bbe'D*Bbe*J*W(k)*W(l)+kem;                 %Bending element stiffness matrix

            %load vector
            Load = zeros(20,1);
            Load(9:20,1):inv(A)*Psi';
            Fe=Load*q0*J*W(k)*W(l);
        end
    end
    Ke=kem+keb;
    for j=1:dofe
        for k=1:dofe
          KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+Ke(j,k);
        end
        FG(CONN(i,j),1)=FG(CONN(i,j),1)+Fe(j,1);
    end
end

%{
%% Preparing matrices for reaction force analysis
FGR = FG;
KGR = KG;

%% Application of Boundary conditions
for kk = [1;2;7;8;13;14]
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end

%% Solving 
UG=KG\FG;
%}



