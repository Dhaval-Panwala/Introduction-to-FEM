clc;
clear ;

Nelem = 5; % no. of element 
Length = 1; % Lehgth of the bar
H = Length/Nelem; % length of element
ApproxOr = 2; % Approximating order
NINT = 6;
ndof = ApproxOr + 1; % Degree of Freedom

for i = 1:Nelem+1
    X(i) = (i - 1)*H; %#ok<*SAGROW>
end

Shape_Function_type = 2; % shape function

%% Shape Function
xi=[-0.9324695142 -0.6612093865 -0.2386191861 0.2386191861 0.6612093865 0.9324695142];
wi=[0.1713244924 0.3607615730 0.4679139346 0.4679139346 0.3607615730 0.1713244924];

if Shape_Function_type == 1 % lagrange
    if ApproxOr == 1 % Linear
        for i = 1:NINT
            N(1,:) = (1 - xi) / 2;
            N(2,:) = (1 + xi) / 2;

            DN(1,:) = -1 / 2 + xi*0;
            DN(2,:) = 1 / 2 + xi*0;
        end

    elseif ApproxOr == 2 % Quadratic
        for i = 1:NINT
            N(1,:) = (xi.*(xi - 1))/2;
            N(2,:) = 1 - xi.^2;
            N(3,:) = (xi.*(xi + 1))/2;

            DN(1,:) = xi - 1 / 2;
            DN(2,:) = -2 * xi;
            DN(3,:) = xi + 1 / 2;
        end

    elseif ApproxOr == 3 % Cubic
        for i = 1:NINT
            N(1,:) = (-(9*xi.^3) + (9*xi.^2) + xi - 1) / 16;
            N(2,:) = ((27*xi.^3) - (9*xi.^2) - (27*xi) + 9) / 16;
            N(3,:) = (-(27*xi.^3) - (9*xi.^2) + (27*xi) + 9) / 16;
            N(4,:) = ((9*xi.^3) + (9*xi.^2) - xi - 1) / 16;

            DN(1,:) = (-(27*xi.^2) + (18*xi) + 1) / 16;
            DN(2,:) = ((81*xi.^2) - (18*xi) - 27) / 16;
            DN(3,:) = (-(81*xi.^2) - (18*xi) + 27) / 16;
            DN(4,:) = ((27*xi.^2) + (18*xi) - 1) / 16;
        end
    elseif ApproxOr == 4 % Quartic
        for i = 1:NINT
            N(1,:) = ((4*xi.^4) - (4*xi.^3) - xi.^2 + xi) / 6;
            N(2,:) = (-(8*xi.^4) + (4*xi.^3) + (8*xi.^2) - (4*xi)) / 3;
            N(3,:) = 4*xi.^4 - 5*xi.^2 + 1;
            N(4,:) = (-(8*xi.^4) - (4*xi.^3) + (8*xi.^2) + (4*xi)) / 3;
            N(5,:) = ((4*xi.^4) + (4*xi.^3) - xi.^2 - xi) / 6;

            DN(1,:) = ((16*xi.^3) - 12*xi.^2 - 2*xi + 1) / 6;
            DN(2,:) = (-(32*xi.^3) + 12*xi.^2 + (16*xi) - 4) / 3;
            DN(3,:) = 16*xi.^3 - 10*xi;
            DN(4,:) = (-(32*xi.^3) - 12*xi.^2 + (16*xi) + 4) / 3;
            DN(5,:) = ((16*xi.^3) + 12*xi.^2 - 2*xi - 1) / 6;
        end
    end
elseif Shape_Function_type == 2 % Heirarchiel
    N(1,:) = (1 - xi) / 2;
    N(2,:) = 3 * (xi.^2 - 1) / (2 * sqrt(6));
    N(3,:) = 5 * (xi.^3 - xi) / (2 * sqrt(10));
    N(4,:) = 7 * (5*xi.^4 - 6*xi.^2 + 1) / (8 * sqrt(14));
    N(5,:) = (1 + xi) / 2;


    N = [N(1:ndof-1,:);N(5,:)];


    DN(1,:) = -1 / 2 + xi*0;
    DN(2,:) = 3 * xi / sqrt(6);
    DN(3,:) = 5 * (3*xi.^2 - 1) / (2 * sqrt(10));
    DN(4,:) = 7 * (5*xi.^3 - 3*xi) / (2 * sqrt(14));
    DN(5,:) = 1 / 2 + xi*0;


    DN = [DN(1:ndof-1,:);DN(5,:)];
end

% Taking Sai and Dsai
for i = 1:NINT
    for j = 1:ApproxOr + 1

        psi(j,i) = N(j,i);
        dpsi(j,i) = DN(j,i);
    end
end

% initiallizing K matrix and F vector  
Ke = zeros(ndof,ndof);
Fe = zeros(1,ndof);

% dof for Global K and F
ntotdof = (Nelem*ApproxOr)+1;

GlobF=zeros(ntotdof,1);
GlobK=zeros(ntotdof,ntotdof);

a = zeros(3);
b = zeros(3);
c = zeros(3);

%% Values of AE,T,and C
a(1) = 1;
b(2) = 1;
c(2) = 1;

for IEL = 1:Nelem
    for i = 1:ndof
        Fe(i) = 0;
            for j = 1:ndof
                Ke(i,j) = 0;
            end
    end

    x1 = X(IEL);
    x2 = X(IEL + 1);
    Jac = (x2 - x1) / 2;

    for l = 1:NINT
        xxi = (x1*((1-xi(l))/2)) + (x2*((1+xi(l))/2));
        AE(l) = a(1) + a(2).*xxi + a(3)*(xxi.^2);
        T(l) = b(1) + b(2).*xxi + b(3)*(xxi.^2);
        %T(l) = sin((pi)*xxi);
        C(l) = c(1) + c(2).*xxi + c(3)*(xxi.^2);

        for i = 1:ndof
            Fe(i) = Fe(i) + T(l) * psi(i,l) * wi(l) * Jac;
            for j = 1:ndof
                Ke(i,j) = Ke(i,j) + ((AE(l) .* dpsi(i,l) .* dpsi(j,l) .* wi(l)) ./ Jac) + (C(l)*psi(i,l)*psi(j,l)*wi(l)*Jac);
            end 
        end
    end

for Irow = 1:ndof

    Irowglo = ieldofs(Irow,IEL,Nelem,ApproxOr,ndof);
    GlobF(Irowglo) = GlobF(Irowglo) + Fe(Irow);
    for Icol = 1:ndof
        Icolglo = ieldofs(Icol,IEL,Nelem,ApproxOr,ndof);
        GlobK(Irowglo,Icolglo) = GlobK(Irowglo,Icolglo) + Ke(Irow,Icol);
    end
end
end

BCvalue_at_LHS = zeros(2,1);
BCvalue_at_RHS = zeros(2,1);

BCatA = 2;
BCvalue_at_LHS(1) = -1;

%% BOundary Conditions at LHS end
if BCatA == 1
    u1 = BCvalue_at_LHS(1);
    for j = 1:ntotdof
        GlobF(j) = GlobF(j)-GlobK(j,1)*BCvalue_at_LHS(1);
    end
    GlobF(1)=BCvalue_at_LHS(1);

    for j = 1:ntotdof
        GlobK(1,j)=0;
        GlobK(j,1)=0;
    end
    GlobK(1,1)=1;

elseif BCatA == 2
    GlobF(1) = GlobF(1)-BCvalue_at_LHS(1);

elseif BCatA == 3
    GlobK(1,1) = GlobK(1,1)+BCvalue_at_LHS(1);
    GlobF(1) = GlobF(1)*(BCvalue_at_LHS(1)*BCvalue_at_LHS(2));

end

BCatB = 1;
BCvalue_at_RHS(1) = 0;
%BCvalue_at_RHS(2) = 0;

%% Boundary Condition at RHS end
if BCatB == 1
    uN = BCvalue_at_RHS(1);
    for j = 1:ntotdof
        GlobF(j) = GlobF(j)-GlobK(j,1)*BCvalue_at_RHS(1);
    end
    GlobF(ntotdof)=BCvalue_at_RHS(1);

    for j = 1:ntotdof
        GlobK(ntotdof,j)=0;
        GlobK(ntotdof,1)=0;
    end
    GlobK(ntotdof,ntotdof)=1;

elseif BCatB == 2
    GlobF(ntotdof) = GlobF(ntotdof)+BCvalue_at_RHS(1);

elseif BCatB == 3
    GlobK(ntotdof,ntotdof) = GlobK(ntotdof,ntotdof)+BCvalue_at_RHS(1);
    GlobF(ntotdof) = GlobF(ntotdof)+(BCvalue_at_RHS(1)*BCvalue_at_RHS(2));

end

%% U vector Calculation
U_Vector = (GlobK)\GlobF;

%% Comparison of exact and FEM solution
%%% Exact Solution
syms ue(x) x xplot
DE = -diff(diff(ue,x))+x*ue == x; % differential eqation
du(x) = diff(ue,x,1); % du/dx;
BC = [du(0)== -1, ue(1) == 0]; % boundary conditions
ue(x) = dsolve(DE,BC);

xvec=0:0.001:1;
ue_sol=subs(ue(x),x,xvec);

du_e(x) = diff(ue,x,1);

x_fem_interp = linspace(0, 1, length(U_Vector));
u_fem_interp = interp1(linspace(0, 1, length(U_Vector)), U_Vector, x_fem_interp);
figure;
plot(xvec,ue_sol, 'b.', x_fem_interp, u_fem_interp, 'ro-');
xlabel('x');
ylabel('u(x)');
legend('Exact Solution', 'FEM Solution');

%% Error Plot
x_fem_interp = linspace(0, 1, length(U_Vector));
u_fem_interp = interp1(linspace(0, 1, length(U_Vector)), U_Vector, x_fem_interp);

% Calculate error
error = ue_sol - interp1(x_fem_interp, u_fem_interp, xvec);

% Plot error
figure;
plot(xvec, error, 'b-', 'LineWidth', 1, 'MarkerSize', 5); 
xlabel('Bar length');
ylabel('Error');
title('Errors ');

%% derivative plot
xvec=0:0.001:1;
du_dx=subs(du_e(x),x,xvec);

% Calculate FEM derivative
x_fem_interp = linspace(0, 1, length(U_Vector));
du_fem_interp = zeros(size(x_fem_interp));
for i = 2:length(x_fem_interp)
    du_fem_interp(i) = (U_Vector(i) - U_Vector(i-1)) / (x_fem_interp(i) - x_fem_interp(i-1));
end

% Plot derivative comparison
figure;
plot(xvec, du_dx, 'b-', x_fem_interp, du_fem_interp, 'o-');
xlabel('Bar length');
ylabel('du/dx');
legend('Exact Derivative', 'FEM Derivative');
title('Derivative Comparison');
%% Function
function IELDOFS = ieldofs(L,K,NELEM,IPVAL,NDOF)
    
    ITEMP = (K -1)*IPVAL;
    JC = ITEMP + L;
    IELDOFS = JC;
end

