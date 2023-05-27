% Assignment 3 - Question 2
clear
% Disclaimer - I have changed the order of the elements and renamed them in
% anti-clockwise order in order to use the regular convention.

%% defining Domain & Grid

points = [0 3 3 0;-1 -1 1 1];

x_div = 21;
y_div = 21;

x_val = linspace(0,3,x_div);
y_val = linspace(-1,1,y_div);

[X,Y] = meshgrid(x_val,y_val);

%% Calculating the interpolated Values

Ae = (points(1,2) - points(1,1))*(points(2,3) - points(2,2)); 

u_inter = zeros(size(X));
v_inter = zeros(size(X));

for i=1:x_div
    for j = 1:y_div

    x = X(j,i);
    y = Y(j,i);
    N1 = (1/Ae) * (x - points(1,2))*(y - points(2,4));
    N2 = -(1/Ae) * (x - points(1,1))*(y - points(2,4));
    N3 = (1/Ae) * (x - points(1,1))*(y - points(2,1));
    N4 = -(1/Ae) * (x - points(1,2))*(y - points(2,1));    
    N = [N1 0 N2 0 N3 0 N4 0;0 N1 0 N2 0 N3 0 N4];
    disp_ele = [-5 5 10 10 15 -10 5 0];
    fval = N*disp_ele';
    u_inter(j,i) = fval(1);
    v_inter(j,i) = fval(2);
    end
end

%% Plotting the Contours

% quiver(X,Y,u_inter,v_inter)
% 
% xlabel('X-Axis')
% ylabel('Y-Axis')
% title('Vector Plot')

% Components of Displacement at the Center of the Element

x_displacement = u_inter(11,11)
y_displacement = v_inter(11,11)

% Total Displacement Calculated Nodewise(Anticlockwise) in tot_disp array.

z_inter = sqrt(u_inter.^2 + v_inter.^2)
