%% Assignment 2
clear

%% defining Domain & Grid

points = [1,3,3,1;1 1 2 2];

x_div = 20;
y_div = 20;

x_val = linspace(1,3,x_div);
y_val = linspace(1,2,y_div);

[X,Y] = meshgrid(x_val,y_val);

%% Calculating Actual Values

ux_act = (3/8)*X - (2/3)*(X./Y.^2);
uy_act = 4./X.^3 - Y/8;

z_act = sqrt(ux_act.^2 + uy_act.^2);

%% Calculating the interpolated Values

Ae = (points(1,2) - points(1,1))*(points(2,3) - points(2,2)); 

u_inter = zeros(size(X));
for i=1:x_div
    for j = 1:y_div
    x = X(j,i);
    y = Y(j,i);
    N1 = (1/Ae) * (x - points(1,2))*(y - points(2,4));
    N2 = -(1/Ae) * (x - points(1,1))*(y - points(2,4));
    N3 = (1/Ae) * (x - points(1,1))*(y - points(2,1));
    N4 = -(1/Ae) * (x - points(1,2))*(y - points(2,1));    
    N = [N1 0 N2 0 N3 0 N4 0;0 N1 0 N2 0 N3 0 N4];
    disp_ele = [ux_act(1,1) uy_act(1,1) ux_act(1,x_div) uy_act(1,x_div) ux_act(y_div,x_div) uy_act(y_div, x_div) ux_act(y_div,1) uy_act(y_div,1)];
    fval = N*disp_ele';
    z_inter = sqrt(fval(1).^2 + fval(2).^2);
    u_inter(j,i) = z_inter;
    end
end

%% Plotting the Contours

hold on;
contour(X,Y,u_inter,25,'-.')
contour(X,Y,z_act,25,'--')
hold off

colorbar()
xlabel('X-Axis')
ylabel('Y-Axis')
title('Both Contour Plot')