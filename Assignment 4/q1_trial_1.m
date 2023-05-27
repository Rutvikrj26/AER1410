% Final Script
% Because of Symmetric arrangement, I am only considering top-right part
clear all; clc
%% Basic Inputs

tot_stress = [];

t = 10; sigma = 120; % [ All Units in mm, Kn, GPa]

tot_disp = [];

% Geometry & Mesh Sizings

lx = 200; ly = 100;
nx = 2; ny = 2;

%% Setting up the Geometry

horizontal_points = linspace(0,lx,nx);
vertical_points = linspace(0,ly,ny);

% Mesh Generation
[X,Y] = meshgrid(horizontal_points,vertical_points);

% Generate Node Coordiantes

% Calculating Mesh Details
total_elements = (nx-1)*(ny-1);
total_nodes = nx*ny;
bottom_nodes = 1:ny;
top_nodes = linspace((ny-1)*nx + 1 ,nx*ny,nx);
side_nodes = nx*[0:(ny-1)] + 1;
% Calculating the Coordinates for all nodes
all_nodes = zeros(total_elements,4);
individual_nodes = [1,2,2+nx,1+nx];
n = 1;
for i = 1:ny-1
for j = 1:nx-1
all_nodes(n,:) = individual_nodes;
element_X_coord(n,:) = [X(i,j),X(i,j+1),X(i+1,j+1),X(i+1,j)];
element_Y_coord(n,:) = [Y(i,j),Y(i,j+1),Y(i+1,j+1),Y(i+1,j)];
individual_nodes = individual_nodes+1;
n = n+1;
end
individual_nodes = individual_nodes+1;
end
%% FEA Part
% Global Stiffness Matrix
Kg = zeros(2*total_nodes);
for i = 1:total_elements
L = gather(total_nodes,all_nodes,i);
Kg = Kg + L'*stiffness(element_X_coord(i,:),element_Y_coord(i,:))*L*t;
end
% Formulaitng the Force Vecor
F = zeros(2*total_nodes,1);
F(2*nx*ny) = 10/sqrt(2);
F(2*nx*ny - 1) = 10/sqrt(2);

% Applying Boundary Conditions

fixeddofs   = sort([2*side_nodes 2*side_nodes-1]);
alldofs     = 1:2*(ny)*(nx);
freedofs    = setdiff(alldofs,fixeddofs);

U(freedofs,:) = Kg(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;

R = Kg*U;
dx = U(1:2:end); dy = U(2:2:end);
tot_def = sqrt(dx.^2 + dy.^2);
loc=1;
for i = 1:ny
    for j =1:nx
        X_mod(i,j) = X(i,j)+ 100*dx(loc); % Applied a Scaling Factor of 100
        Y_mod(i,j) = Y(i,j)+ 100*dy(loc); % Applied a Scaling Factor of 100
        loc = loc+1;
    end
end

Sg_elements = zeros(total_elements,2);
Sg_nodes = zeros(total_nodes,2);
node_counter = zeros(total_nodes,1);
for i = 1:total_elements
element_nodes = all_nodes(i,:);
Q = zeros(8,1);
for nodes =1:4
Q(2*nodes-1) = dx(element_nodes(nodes));
Q(2*nodes) = dy(element_nodes(nodes));
end
Sg_elements(i,:) = stress(element_X_coord(i,:),element_Y_coord(i,:),Q);
for j = 1:4
curr_node = all_nodes(i,j);
Sg_nodes(curr_node,1) = Sg_nodes(curr_node,1) + Sg_elements(i,1);
Sg_nodes(curr_node,2) = Sg_nodes(curr_node,2) + Sg_elements(i,2);
node_counter(curr_node) = node_counter(curr_node) + 1;
end
end
Sg_nodes = Sg_nodes./node_counter;
max(Sg_nodes(:,2))
%% Contours and Plots

figure(1)
hold on
 for i = 1:ny
 plot(X(i,:),Y(i,:),'k')
 plot(X_mod(i,:),Y_mod(i,:),'b')
 end
plot(X,Y,'k')
plot(X_mod,Y_mod,'b')

title('Mesh with & w/o Deflection(100x)')
hold off

von_misses = transpose(reshape(Sg_nodes(:,1),ny,nx));
figure(4)
contourf(X,Y,von_misses,20)

title('Von-Misses Stress')
colorbar()

tot_stress = [tot_stress (Sg_nodes(1,2))];

%% Functions for Stiffness & Gather Matrix Calculations
function K = stiffness(X,Y)
E = 120; % [in GPa]
nu = 0.25;
D = (E/((1+nu)*(1-2*nu)))*[1-nu nu 0; nu 1-nu 0;0 0 0.5-nu];
coord = [X',Y'];
K = zeros(8,8);
for i = 1:2
for j = 1:2
eta = (2*i-3)/sqrt(3);
zeta = (2*j-3)/sqrt(3);
J = (1/4)*[eta-1 1-eta 1+eta -eta-1; zeta-1 -zeta-1 1+zeta 1-zeta]*coord;
H = (1/4)*[eta-1 1-eta 1+eta -eta-1; zeta-1 -zeta-1 1+zeta 1-zeta];
H = J\H;
H = [H(1,1) 0 H(1,2) 0 H(1,3) 0 H(1,4) 0; 0 H(2,1) 0 H(2,2) 0 H(2,3) 0 H(2,4); H(2,1) H(1,1) H(2,2) H(1,2) H(2,3) H(1,3) H(2,4) H(1,4)];
K = K + det(J)*H'*D*H;
end
end
end

function S = stress(X,Y,Q)
E = 120; % [GPa]
nu = 0.25;
D = (E/((1+nu)*(1-2*nu)))*[1-nu nu 0; nu 1-nu 0;0 0 0.5-nu];
coord = [X',Y'];
stress_temp = zeros(3,1);
for i = 1:2
for j = 1:2
eta = (2*i-3)/sqrt(3);
zeta = (2*j-3)/sqrt(3);
J = (1/4)*[eta-1 1-eta 1+eta -eta-1; zeta-1 -zeta-1 1+zeta 1-zeta]*coord;
H = (1/4)*[eta-1 1-eta 1+eta -eta-1; zeta-1 -zeta-1 1+zeta 1-zeta];
H = J\H;
H = [H(1,1) 0 H(1,2) 0 H(1,3) 0 H(1,4) 0; 0 H(2,1) 0 H(2,2) 0 H(2,3) 0 H(2,4); H(2,1) H(1,1) H(2,2) H(1,2) H(2,3) H(1,3) H(2,4) H(1,4)];
stress_temp = stress_temp + 0.25*D*H*Q; % Sigma XX, Sigma YY, Sigma XY
end
end
stress_xx = stress_temp(1);
stress_yy = stress_temp(2);
stress_xy = stress_temp(3);
stress_vm = sqrt(stress_xx^2 + stress_yy^2 + 3*stress_xy^2 - stress_xx*stress_yy);
S = [stress_vm stress_yy];
end

function L = gather(total_nodes,elements,n)
L = zeros(8,2*total_nodes);
L(1:2,2*elements(n,1)-1:2*elements(n,1)) = eye(2);
L(3:4,2*elements(n,2)-1:2*elements(n,2)) = eye(2);
L(5:6,2*elements(n,3)-1:2*elements(n,3)) = eye(2);
L(7:8,2*elements(n,4)-1:2*elements(n,4)) = eye(2);
end