%% Vibrational Analysis of plate using FSDT quadrilateral elements
clc                                            %clears screeen
clear all                                      %clears workspace

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
G = 90E9;                                       %Shear Modulus
nu = 0.25;                                      %Poisson's ratio
Rho = 7.8;                                      %Density
t = 0.1;                                        %thickness of plate
lx = 0.5;                                         %length of element in x-direction
ly = 0.5;                                         %ength of element in y-direction
%q0 = 100;                                       %uniformly distributed load
Nxx=0;                                          %Load in x-direction
Nyy=0;                                          %Load in y-direction
Nxy=0;                                          %Load in z-direction

%% Assembly of element stiffness matrices
% 4 point formula
Gp=[-0.8611 -0.3399 0.3399 0.8611];
Gp2 = [-0.57735,0.57735];                       %Gauss Points
W=[0.3478 0.652 0.652 0.3478];
W2 = [1,1];

p=size(Gp,2);                                   %No. of columns in Gp
q=size(Gp,2);
o=size(Gp2,2);
m=size(Gp2,2);
D1=(E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2];
D2 = 5*[G,0;0,G]/6;
KG=zeros(tdof,tdof);
MG = zeros(tdof);
for i=1:ne                                      %Membrane stiffness
    de=[coord(NCONN(i,1),1),coord(NCONN(i,1),2);coord(NCONN(i,2),1),coord(NCONN(i,2),2);coord(NCONN(i,3),1),coord(NCONN(i,3),2);coord(NCONN(i,4),1),coord(NCONN(i,4),2)];
    kme=zeros(dofe,dofe);
    kbe=zeros(dofe,dofe);
    Kshe=zeros(dofe);
    Ke=zeros(dofe);
    Me = zeros(dofe);
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

            %Membrane bending matrix
            Bme = zeros(3,dofe);
            Bbe = zeros(3,dofe);
            for a=1:nne
                Bme(1,5*a-4) = phi_x(1,a);
                Bme(2,5*a-3) = phi_y(1,a);
                Bme(3,5*a-3) = phi_x(1,a);
                Bme(3,5*a-4) = phi_y(1,a);
                Bbe(1,5*a-1) = phi_x(1,a);
                Bbe(2,5*a) = phi_y(1,a);
                Bbe(3,5*a) = phi_x(1,a);
                Bbe(3,5*a-1) = phi_y(1,a);
            end
            kme=t*Bme'*D1*Bme*J*W(k)*W(l)+kme;
            kbe = t^3/12*Bbe'*D1*Bbe*J*W(k)*W(l)+kbe;

            %Mass Matrix
            S=zeros(5,20);
            Y=zeros(5,5);
            for ii=1:5
                S(ii,4*(ii-1)+1:4*(ii-1)+4)=phi;
            end
            I0=Rho*t;
            I1=0;
            I2=Rho*t^3/12;

            P=[I0,0,0,I1,0;0,I0,0,0,I1;0,0,I0,0,0;I1,0,0,I2,0;0,I1,0,0,I2];

            Me=S'*P*S*J*W(k)*W(l)+Me;

            %{
            %load Vector
            psi=zeros(1,dofe);
            for a=1:nne
                psi(1,5*a-2) = phi(1,a);
            end

            Fe=psi'*q0*J*W(k)*W(l);
            %}
            end
        end
        
    %start second loop for shear
    for k = 1:o
        for l = 1:m
            e=Gp2(k);    % Taking Gauss points into e and n
            n=Gp2(l);
    
            phi = 0.25*[(1-e)*(1-n), (1+e)*(1-n), (1+e)*(1+n), (1-e)*(1+n)];
            phi_z = 0.25*[-(1-n), (1-n), (1+n), -(1+n)];
            phi_n = 0.25*[-(1-e), -(1+e), (1+e), (1-e)];
            X = [de(1,1); de(2,1); de(3,1); de(4,1)];
            Y = [de(1,2); de(2,2); de(3,2); de(4,2)];
    
 
            %Assembling phi 
            x_z=0.25*[-(1-n), (1-n), (1+n), -(1+n)]*X;   
            x_n=0.25*[-(1-e), -(1+e), (1+e), (1-e)]*X;
            y_z=0.25*[-(1-n), (1-n), (1+n), -(1+n)]*Y;
            y_n=0.25*[-(1-e), -(1+e), (1+e), (1-e)]*Y;
                
            J=x_z*y_n-x_n*y_z;
            phi_x = y_n/J*phi_z-y_z/J*phi_n;
            phi_y = -x_n/J*phi_z+x_z/J*phi_n;
            Bshe = zeros(2,dofe);
            
            for a=1:nne
                Bshe(1,5*a-2)=phi_x(1,a);
                Bshe(1,5*a-1)=phi(1,a);
                Bshe(2,5*a-2)=phi_y(1,a);
                Bshe(2,5*a)=phi(1,a);
            end
            Kshe=t*Bshe'*D2*Bshe*J*W(k)*W(l);

        end
    end
    Ke=kme+kbe+Kshe;
    for j=1:dofe
        for k=1:dofe
          KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+Ke(j,k);
          MG(CONN(i,j), CONN(i,k)) = MG(CONN(i,j), CONN(i,k)) + Me(j,k);
        end
    end
end

%% Preparing matrices for reaction force analysis
KGR = KG;
MGR = MG; 
UGR = zeros(tdof,1);

%% Application of Boundary conditions
c=1;                                                                %element number count
cc=1;                                                               %Boundary element index
boundary=zeros(20+6*(x-1)+6*(y-1),1);
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
    MG(kk,:) = [];
    MG(:,kk) = [];
end


%% Solving
[eigvect, eigval] = eig(KG, MG);                                        %Eigenvalue solving
[omega, index] = sort(sqrt(diag(eigval)));
freq = omega/(2*pi);

%% Post Processing
UG = KG-omega(1)^2*MG;
UGR = zeros(tdof);

deflection=zeros(size(omega,1),size(omega,1));
for i=1:1:size(omega,1)
    deflection(:,i)=eigvect(:,i)/eigvect(size(omega,1),1);
end

count1=1;
count2=1;
for jj=1:tdof
    if jj==boundary(count1)
        count1=count1+1;
    else
        UGR(jj)=deflection(3,count2);
        count2=count2+1;
    end
end

ca=1;
zcoord=zeros(x+1,y+1);
for a=1:y+1
    for b=1:x+1
        zcoord(a,b)=UGR(5*ca-2);
        ca=ca+1;
    end
end

%% Output
disp('The 2 lowest natural frequencies are')
disp(freq(1:2,1))