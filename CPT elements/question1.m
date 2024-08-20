%% Static Analysis of plate using CPT quadrilateral elements
clc                             %clears screeen
clear all                       %clears workspace

%% Meshing
x=2;                                            %Number f nodes along x-axis
y=2;                                            %Number of nodes along y-axis
nne = 4;                                        %Numeber of nodes per element
dofn = 5;                                       %Degrees of freedom per node
ne = x*y;                                       %Number of elements
nn = (x+1)*(y+1);                               %Number of nodes
dofe = nne*dofn;                                %Degrees of freedom per element
tdof = dofn*nn;                                 %Total Degrees of Freedom
NCONN=zeros(ne,nne);                            %Declarin NCONN
count=1;
for a=1:y
    for i=1:x
        if a==1
            NCONN(count,:)= [i,i+1,x+i+2,x+i+1];
        else
            NCONN(count,:)= [x*a+i-1,x*a+i,x*a+i-1+x+2,x*a+i+x+1];
        end
        count=count+1;
    end
end
CONN=zeros(ne,dofe);                            %Connectivity Matrix
for i=1:ne
    for j=1:nne
        for k=1:dofn
            CONN(i,(j-1)*5+k) = (NCONN(i,j)-1)*5+k;
        end
    end
end
coord = zeros(nn,2);                             %declaring coord matrix
count2=1;
for i=1:y+1
    for j=1:x+1
        coord(count2,:)=[(i-1),(j-1)];
        count2=count2+1;
    end
end

%% Geometry and Material Properties
Ae = 0.01;                                      %Area of Element
E = 200E09;                                     %Young's modulus
nu = 0.3;                                      %Poisson's ratio
t = 0.05;                                       %thickness of plate
lx = 1;                                         %length of element in x-direction
ly = 1;                                         %ength of element in y-direction
q0 = 100;                                 %uniformly distributed load

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
    de=[coord(NCONN(i,1),1),coord(NCONN(i,1),2);coord(NCONN(i,2),1),coord(NCONN(i,2),2);coord(NCONN(i,3),1),coord(NCONN(i,3),2);coord(NCONN(i,4),1),coord(NCONN(i,4),2)];
    kme=zeros(dofe,dofe);
    kbe=zeros(dofe,dofe);
    Kb=zeros(dofe);
    Fe = zeros(dofe,1);
    for k = 1:p
        for l = 1:q
            e=Gp(k);    % Taking Gauss points into e and n
            n=Gp(l);

            phi = 0.25*[(1-e)*(1-n), (1+e)*(1-n), (1+e)*(1+n), (1-e)*(1+n)];
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
                Bme(1,5*a-4) = phi_x(1,a);
                Bme(2,5*a-3) = phi_y(1,a);
                Bme(3,5*a-3) = phi_x(1,a);
                Bme(3,5*a-4) = phi_y(1,a);
            end
            kme =t*Bme'*D*Bme*J*W(k)*W(l)+kme;                      %Membrane element stiffness matrix

            %bending stiffness matrix
            A = [1,-1,-1,1,1,1,-1,-1,-1,-1,1,1; 1,1,-1,1,-1,1,1,-1,1,-1,-1,-1; 1,1,1,1,1,1,1,1,1,1,1,1; 1,-1,1,1,-1,1,-1,1,-1,1,-1,-1; 0,1/lx,0,-2/lx,-1/lx,0,3/lx,2/lx,1/lx,0,-3/lx,-1/lx; 0,1/lx,0,2/lx,-1/lx,0,3/lx,-2/lx,1/lx,0,-3/lx,-1/lx; 0,1/lx,0,2/lx,1/lx,0,3/lx,2/lx,1/lx,0,3/lx,1/lx; 0,1/lx,0,-2/lx,1/lx,0,3/lx,-2/lx,1/lx,0,3/lx,1/lx; 0,0,1/ly,0,-1/ly,-2/ly,0,1/ly,2/ly,3/ly,-1/ly,-3/ly; 0,0,1/ly,0,1/ly,-2/ly,0,1/ly,-2/ly,3/ly,1/ly,3/ly; 0,0,1/ly,0, 1/ly,2/ly,0,1/ly,2/ly,3/ly,1/ly,3/ly; 0,0,1/ly,0,-1/ly,2/ly,0,1/ly,-2/ly,3/ly,-1/ly,-3/ly];
            Psi = [1,e,n,e^2,e*n,n^2,e^3,e^2*n,e*n^2,n^3,e^3*n,e*n^3];
            Psi_z = [0,1,0,2*e,n,0,3*e^2,2*e*n,n^2,0,3*e^2*n,n^3];
            Psi_n = [0,0,1,0,e,2*n,0,2*e^2,2*e*n,3*n^2,e^3,3*n^2*e];
            Psi_zz = [0,0,0,2,0,0,6*e,2*n,0,0,6*e*n,0];
            Psi_nn = [0,0,0,0,0,2,0,0,2*e,6*n,0,6*n*e];
            Psi_nz = [0,0,0,0,1,0,0,2*e,2*n,0,3*e^2,3*n^2];
            
            B = Psi/A;                                         %defining Beta matrices
            Jw = [x_z,y_z,0,0,0; x_n,y_n,0,0,0; 0,0,x_z^2,y_z^2,2*x_z*y_z; 0,0,x_n^2,y_n^2,2*x_n*y_n; 0,0,x_n*x_z,y_n*y_z,x_z*y_n+x_n*y_z];
            Psicart = Jw\[Psi_z;Psi_n;Psi_zz;Psi_nn;Psi_nz];
            Bxx = Psicart(3,1:12)/A;
            Byy = Psicart(4,1:12)/A;
            Bxy = Psicart(5,1:12)/A;
            
            Bbe = zeros(3,20);                                  
            for a=1:nne
                Bbe(1,5*a-2) = Bxx(1,3*a-2);
                Bme(1,5*a-1) = Bxx(1,3*a-1);
                Bme(1,5*a) = Bxx(1,3*a);
                Bbe(2,5*a-2) = Byy(1,3*a-2);
                Bme(2,5*a-1) = Byy(1,3*a-1);
                Bme(2,5*a) = Byy(1,3*a);
                Bbe(3,5*a-2) = 2*Bxy(1,3*a-2);
                Bme(3,5*a-1) = 2*Bxy(1,3*a-1);
                Bme(3,5*a) = 2*Bxy(1,3*a);
            end
            
            kbe = t^3/12*Bbe'*D*Bbe*J*W(k)*W(l)+kbe;                 %Bending element stiffness matrix

            %load vector
            Load = zeros(dofe,1);
            Loading= A\Psi';
            for a=1:nne
                Load(5*a-2,1) = Loading(3*a-2,1);
                Load(5*a-1,1) = Loading(3*a-1,1);
                Load(5*a,1) = Loading(3*a,1);
            end
            Fe=Load*q0*J*W(k)*W(l);
        end
    end
    Ke=kme+kbe;
    % Fe = t*le/2*[0; 0; Tx(i); 0; Tx(i); 0; 0; 0]; 
    for j=1:dofe
        for k=1:dofe
          KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+Ke(j,k);
        end
        FG(CONN(i,j),1)=FG(CONN(i,j),1)+Fe(j,1);
    end
end

%% Preparing matrices for reaction force analysis
FGR = FG;
KGR = KG;


%% Application of Boundary conditions
c=1;                                                                %element number count
cc=1;                                                               %Boundary element index
boundary=zeros(16+6*(x-1)+6*(y-1),1);
for i=1:y+1
    for j=1:x+1
        if i==1 && j==1 
            boundary(cc:cc+4,1)=[5*c-4;5*c-3;5*c-2;5*c-1;5*c];
            cc=cc+5;
        elseif i==y+1 && j==1
            boundary(cc:cc+4,1)=[5*c-4;5*c-3;5*c-2;5*c-1;5*c];
            cc=cc+5;
        elseif i==1 && j==x+1
            boundary(cc:cc+4,1)=[5*c-4;5*c-3;5*c-2;5*c-1;5*c];
            cc=cc+5;
        elseif i==y+1  && j==x+1
            boundary(cc:cc+4,1)=[5*c-4;5*c-3;5*c-2;5*c-1;5*c];
            cc=cc+5;
        elseif i==1 || i==y+1
            boundary(cc:cc+2,1)=[5*c-4;5*c-2;5*c-1];
            cc=cc+3;
        elseif j==1 || j==x+1
            boundary(cc:cc+2,1)=[5*c-3;5*c-2;5*c];
            cc=cc+3;
        end
        c=c+1;
    end
end
for kk=boundary
    KG(kk,:) = [];
    KG(:,kk) = [];
    FG(kk,:) = [];
end


%% Solving 
UG=KG\FG

%{
%% Post Processing
UGD = zeros(nn,1); 
UGR = zeros(tdof,1);
UGR(3:6,1) = UG(1:4,1);
UGR(9:12,1) = UG(5:8,1);
UGR(15:18,1) = UG(9:12,1);
for p = 1:nn
    UGD(p,1) = ((UGR(2*p-1,1))^2+(UGR(2*p,1))^2)^(1/2);
end
RF = KGR*UGR-FGR;
%}

%% Output

