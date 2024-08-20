%% Non-linear analysis of a straight bar
clc                                     %Clears command window
clear all                               %Clears workspace

%% Geometry and material properties
L = 1;                                  %Length of bar
E = 100E09;                             %Young's Modulus
nu = 0.25;
A = 0.0002;                             %Area of Cross section
t = 0.1;                                %Thickness
P = 100;                                %concentrated load
q0 = 0;                                 %uniform load

%% Meshing
ne = 1;                                 %Number of elements
nne = 2;                                %Number of nodes per element
nn = ne+1;                              %Number of nodes
dofn = 1;                               %Degrees of freedom per node
dofe = dofn*nne;                        %Degrees of freedom per element
tdof = dofn*nn;                         %Total degrees of freedom
le = L/ne;                              %length of each element
CONN = zeros(ne, dofe);                 %Initializing connectivity matrix
for ii=1:ne
    %CONN(ii,:) = ii:ii+1;
    CONN(ii,:) = (dofe-1)*(ii-1)+1: dofe + (dofe-1)*(ii-1);
end

%% Initial guess for loop
delu = [0;0];
sumsq1 = 0;
sumsq2 = 0;
displacement = zeros(11,1);
DELF = zeros(11,1);
z=0;

for UI = 0:0.0005:0.005
    conv=10;
    ui=[0;UI];
    while conv>0.001
        KG = zeros(tdof);                                               %Initializing global stiffness matrix
        FG = zeros(tdof,1);                                             %Initializing global load matrix   
        KTG = zeros(tdof);                                              %Inititalizing tangent matrix\
            
        D=E;
            
        GP=[-0.8611 -0.3399 0.3399 0.8611];
        W=[0.3478 0.652 0.652 0.3478];
            
        p = size(GP,1);
        q = size(GP,1);
            
        for i = 1:ne
            Ke = zeros(dofe);
            Kt = zeros(dofe);
            Fe = zeros(dofe,1);
            Kg = zeros(dofe);
            
            phi_z = 1/2*[-1,1];
            Bl = 2/le*phi_z;
            Bnl = (ui)'*4/le^2*phi_z'*phi_z;
                    
            Kle = (Bl')*D*Bl;
            Knle1 = (Bnl')*D*Bl;
            Knle2 = (Bl')*D*Bnl;
            Knle3 = (Bnl')*D*Bnl;
    
            Nxx = E*(Bl+Bnl/2)*(ui)*le;
            
            Kg = phi_z'*4*Nxx/le^2*phi_z*A +Kg;
            
            Ke = Kle+Knle1+Knle2+Knle3;
            Kt = Kle+Knle1+Knle2+Knle3;%+Kg;
            
            
            for j = 1:dofe
                for k = 1:dofe
                    KG(CONN(i,j), CONN(i,k)) = KG(CONN(i,j), CONN(i,k)) + Ke(j,k);
                    KTG(CONN(i,j), CONN(i,k)) = KTG(CONN(i,j), CONN(i,k)) + Kt(j,k);
                end
                    %FGU(CONN(i,j),1) = FGU(CONN(i,j),1) + FE(j,1);
            end
        end
    
        %% Preaparing solving matrices and applying boundary conditions
        KGR = KG;
        FGR = FG;
        UGR = zeros(tdof,1);
        KG(1,:) = [];
        KG(:,1) = [];
        FG(2,1) = P;
        FG(1,:) = [];
        KTG(1,:) = [];
        KTG(:,1) = [];
    
    
        %% Netwon-Raphson
        DelF = FG-KG*ui(2);
        delu(2) = KTG\DelF;
        ui = ui+delu;
    
        sumsq1 = sumsq1+delu(2)^2;
        sumsq1 = sumsq1+ui(2)^2;
    
        conv = norm(delu);
    
    end
    %disp('next displacement step');
    z=z+1;
    displacement(z) = ui(2);
    DELF(z)=DelF;
end

F=[100,100,100,100,100,100,100,100,100,100,100];

plot(DELF,displacement)

disp('The displacement at the tip is')
disp(displacement(11))
disp('The strain at the tip is')
disp(100/0.0002/E)
disp('The stress at the tip is')
disp(100/0.0002)




