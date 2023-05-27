% Assignment 3 - Question 3

points = [0 3 3 0;-1 -1 1 1];

x_e = points(1,:);
y_e = points(2,:);

de = [-5 5 10 10 15 -10 5 0]';

E = 110e+9;
nu = 0.3;
x = 1;
y = 0.5;

Ae = 

D = (E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2];

H = (1/Ae)*[(y - y_e(4)), 0, -(y - y_e(4)), 0 (y - y_e(1)), 0 , -(y - y_e(1)), 0;
            0, (x - x_e(2)), 0, -(x - x_e(1)), 0, (x - x_e(1)), 0, -(x - x_e(2));
            (x - x_e(2)) (y - y_e(4)) -(x - x_e(1)) -(y - y_e(4)) (x - x_e(1)) (y - y_e(1)) -(x - x_e(2)) -(y - y_e(1))];

strain = H*de
stress = D*strain
