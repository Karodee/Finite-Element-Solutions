%% Newton-Raphson method for non-linear analysis 
clc                                                     %Clears screen
clear all                                               %Clears workspace

%% Defining functions and variables
syms x y;                                               %Symbolic variables
function1 = x^2+y^2-3;                                         %First function
function2 = x*y-1;                                             %Second function
sumsq = 0;                                              %Sum of squares
conv = 4/5;                                 %Conversion criteria expression
count=0;                                                %Jugaad

%% Defining Matrices
m=2; 
n=1; 
F1 = [0;0];                                             %Initial F matrix
delF = zeros(2,1);                                      %Delta F matrix
X = zeros(2,1);                                         %X matrix

%% Iteration loop
for i=1:3
    A = subs(diff(function1,x),{x,y},{m,n});
    B = subs(diff(function1,y),{x,y},{m,n});
    C = subs(diff(function2,x),{x,y},{m,n});
    D = subs(diff(function2,y),{x,y},{m,n});
    f1 = subs(function1,{x,y},{m,n});
    f2 = subs(function2,{x,y},{m,n});

    F1

    Kt = [A,B;C,D]                                         %Tangent matrix
    F = F1-[f1;f2]                                         %Function matrix
    X = Kt\F;

    sumsq = sumsq + (X(1)-m)^2;

    m=m+X(1);
    n=n+X(2);
    count+count+1;
    conv = sumsq/(1+sumsq);
    F1=F;
    disp([m;n])
end

disp('The intersection point is ')
disp([m;n])
    

