clc;
clear;
% Input parameters
Nelement = input('Enter number of elements:'); % number of element 
Nnodes = Nelement + 1; % number of nodes 
ndof = Nnodes*2; % Degree of freedom
E = 2e7; % Young's Modulus( in N/cm2)
H = 1; % Height of beam (in cm)
width = 1; % (in cm)
Ltotal = 10; % length of beam (in cm)
I = (width*H^3)/12; % Moment of inertia (in cm^4)
Lelement = Ltotal/Nelement; % element length

% Ask user for Load type
fprintf('\nSelect Load type:\n');
fprintf('1. Distributed load over beam\n');
fprintf('2. Point load at half\n');
fprintf('3. Point load at end\n');
fprintf('4. Point moment at half\n');
Load = input('Enter Load type (1, 2, 3, or 4): ');

% Ask user for Boundary condition at left
fprintf('\nSelect Boundary condition at left:\n');
fprintf('1. Displacement as a BC\n');
fprintf('2. Slope as a BC\n');
fprintf('3. Slope and Displacement as a BC\n');
fprintf('4. Slope and shear force as Boundary condition\n');
fprintf('5. Displacement and Moment as a BC\n');
fprintf('6. Shear and Moment as BC\n');
fprintf('7. Shear as a BC\n');
fprintf('8. Moment as a BC\n');
BC_left = input('Enter Boundary condition at left (1-8): ');

% Ask user for Boundary condition at right
fprintf('\nSelect Boundary condition at right:\n');
fprintf('1. Displacement as a BC\n');
fprintf('2. Slope as a BC\n');
fprintf('3. Slope and Displacement as a BC\n');
fprintf('4. Slope and shear force as Boundary condition\n');
fprintf('5. Displacement and Moment as a BC\n');
fprintf('6. Shear and Moment as BC\n');
fprintf('7. Shear as a BC\n');
fprintf('8. Moment as a BC\n');
fprintf('9. Free Boundary\n');
BC_right = input('Enter Boundary condition at right (1-9): ');

% Initialize matrices and vectors
Q = zeros(ndof,1);
GlobK = zeros(ndof,ndof);
GlobF = zeros(ndof,1);

% Calculate element stiffness matrix
kelement =((E*I)/Lelement^3)*[ 12, -6*Lelement, -12, -6*Lelement;
            -6*Lelement, 4*Lelement^2, 6*Lelement, 2*Lelement^2;
            -12, 6*Lelement, 12, 6*Lelement;
            -6*Lelement, 2*Lelement^2, 6*Lelement, 4*Lelement^2];

% Assemble stiffness matrix
for i = 1:Nelement
    GlobK(2*i-1 : 2*(i+1), 2*i-1 : 2*(i+1)) = GlobK(2*i-1 : 2*(i+1), 2*i-1 : 2*(i+1)) + kelement;
end

% Calculate element force vector
if Load == 1
    q = input('\nUniformly distributed load value (in N):'); % Distributed load over beam
    qe = q; % distributed load over element
    felement = ((qe*Lelement)/12)*[6; -Lelement; 6; Lelement];
elseif Load == 2
    Q(ndof/2) = input('\nPoint load value at center of the beam length (in N):');% point load at half
    felement = zeros(4,1);
elseif Load == 3
    Q(ndof-1) = input('\nPoint load value at end of the beam length (in N):'); %point load at end
    felement = zeros(4,1);
elseif Load == 4
    Q(ndof/2 + 1) = input('\nPoint Moment value at center of the beam length (in N):');% point moment at half
    felement = zeros(4,1);
end

% Assemble force vector
for i = 1:Nelement
    GlobF(2*i-1) = GlobF(2*i-1) + felement(1); % horizontal displacement of node i
    GlobF(2*i) = GlobF(2*i) + felement(2); % vertical displacement of node i
    GlobF(2*(i+1)-1) = GlobF(2*(i+1)-1) + felement(3); % horizontal displacement of node i+1
    GlobF(2*(i+1)) = GlobF(2*(i+1)) + felement(4); % vertical displacement of node i+1
end

% Boundary Condition at left end side
if BC_left == 1 % Displacement as a BC
    GlobK(1,:) = 0;
    GlobK(:,1) = 0;
    GlobK(1,1) = 1;
    GlobF(1,1) = 0;
elseif BC_left == 2 % Slope as a BC
    GlobK(2,:) = 0;
    GlobK(:,2) = 0;
    GlobF(2,2) = 1;
    GlobF(2,1) = 0;
elseif BC_left == 3 % Slope and Displacement as a BC 
    GlobK(1,:) = 0;
    GlobK(:,1) = 0;
    GlobK(1,1) = 1;
    GlobK(2,:) = 0;
    GlobK(:,2) = 0;
    GlobK(2,2) = 1;
    GlobF(1,1) = 0;
    GlobF(2,1) = 0;
elseif BC_left == 4 % Slope and shear force as Boundary condition;
    S_left = input('\nEnter Shear force at left end side:');
    Q(1,1) = S_left;
    GlobK(2,:) = 0;
    GlobK(:,2) = 0;
    GlobF(2,2) = 1;
    GlobF(2,1) = 0;
elseif BC_left == 5 % displacement and Moment as a BC
    M_left = input('\nEnter Moment at left end side:');
    Q(2,1) = M_left;
    GlobK(1,:) = 0;
    GlobK(:,1) = 0;
    GlobK(1,1) = 1;   
    GlobF(1,1) = 0;
elseif BC_left == 6 % Shear and Moment as BC
    S_left = input('\nEnter Shear force at left end side:');
    Q(1,1) = S_left;
    M_left = input('\nEnter Moment at left end side:');
    Q(2,1) = M_left;
elseif BC_left == 7 % shear as a BC
    S_left = input('\nEnter Shear force at left end side:');
    Q(1,1) = S_left;
elseif BC_left == 8 % Moment as a BC  
    M_left = input('\nEnter Moment at left end side:');
    Q(2,1) = M_left;
end

% Boundary condition at Right side
if BC_right == 1 % Displacement as a BC
    GlobK(ndof-1,:) = 0;
    GlobK(:,ndof-1) = 0;
    GlobK(ndof-1,ndof-1) = 1;
    GlobF(ndof-1 ,1) = 0;
elseif BC_right == 2 % Slope as a BC
    GlobK(ndof,:) = 0;
    GlobK(:,ndof) = 0;
    GlobF(ndof,ndof) = 1;
    GlobF(ndof,1) = 0;
elseif BC_right == 3 % Slope and Displacement as a BC 
    GlobK(ndof-1,:) = 0;
    GlobK(:,ndof-1) = 0;
    GlobK(ndof-1,ndof-1) = 1;
    GlobK(ndof,:) = 0;
    GlobK(:,ndof) = 0;
    GlobK(ndof,ndof) = 1;
    GlobF(ndof-1,1) = 0;
    GlobF(ndof,1) = 0;
elseif BC_right == 4 % Slope and shear force as Boundary condition
    S_right = input('\nEnter Shear Force at Right end side:');
    Q(ndof-1,1) = S_right;
    GlobK(ndof,:) = 0;
    GlobK(:,ndof) = 0;
    GlobF(ndof,ndof) = 1;
    GlobF(ndof,1) = 0;
elseif BC_right == 5 % displacement and Moment as a BC
    M_right = input('\nEnter Moment at left end side:');
    Q(ndof,1) = M_right;
    GlobK(ndof-1,:) = 0;
    GlobK(:,ndof-1) = 0;
    GlobK(ndof-1,ndof-1) = 1;   
    GlobF(ndof-1,1) = 0;    
elseif BC_right == 6 % Shear and Moment as BC
    S_right = input('\nEnter Shear Force at Right end side:');
    Q(ndof-1,1) = S_right;
    M_right = input('\nEnter Moment at Right end side:');
    Q(ndof,1) = M_right;
elseif BC_right == 7 % shear as a BC
    S_right = input('\nEnter Shear Force at Right end side:');
    Q(ndof-1,1) = S_right;
elseif BC_right == 8 % Moment as a BC   
    M_right = input('\nEnter Moment at Right end side:');
    Q(ndof,1) = M_right;
elseif BC_right == 9 % Free Boundary
    R = zeros(ndof,ndof);
    D = zeros(ndof,1);
    GlobK = GlobK + R;
    GlobF = GlobF + D;
end
% Post processing
Force = GlobF + Q;
GlobK = GlobK + eye(size(GlobK))*1e-6;

% Solve the system using mldivide
U = inv(GlobK)*Force;

% Extracting displacement and slope from U
dispplacement = -U(1:2:end);
slope = U(2:2:end);

% Creating a vector of x-coordinates
x = linspace(0, Ltotal, Nnodes);

% Plotting displacement
figure;
plot(x, dispplacement);
hold on
plot(x, slope);
xlabel('Distance from left end (cm)');
ylabel('Slope (rad) and displacement (cm)');
title('Slope and displacement along the Beam');
legend('displacement', 'slope')
hold off
grid on

% Calculate derivatives using finite differences
dx = Ltotal / (Nnodes-1);
dw_dx = zeros(Nnodes-1, 1);
d2w_dx2 = zeros(Nnodes-2, 1);

for i = 1:Nnodes-1
    dw_dx(i) = (U(2*i+1) - U(2*i-1)) / (dx);
end

for i = 1:Nnodes-2
    d2w_dx2(i) = (dw_dx(i+1) - dw_dx(i)) / dx;
end

% Calculate bending moment
M = -E * I * d2w_dx2;


% Calculate shear force
V = zeros(Nnodes-3, 1);
for i = 1:Nnodes-3
    V(i) = -(M(i+1) - M(i)) / dx;
end


% Plot bending moment and shear force
figure;
hold on
plot(x(2:end-2), V, '-b'); % Modify the indexing here
xlabel('Distance from left end (cm)');
ylabel('Bending Moment (Ncm) and Shear Foece (N)');
hold on
plot(x(2:end-1), M, '-r');

legend('Shear Force','Bending Moment');

grid on;hold off;
figure;
stress= H/(2*I)*M;
plot(x(2:end-1), stress, 'b');
title('Bending Stress at top most layer');
ylabel('Bending Stress');
xlabel('Distance from left end (cm)');
grid on;
ylabel('Deflection');