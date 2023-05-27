%% Assignment-3 Q1
clear
% Disclaimer - I have changed the order of the elements and renamed them in
% anti-clockwise order in order to use the regular convention.

E = 110e+9;
nu = 0.3;
points = [0 3 3 0;-1 -1 1 1];

x_e = points(1,:);
y_e = points(2,:);
% Plane Stress Conditions
D = (E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2];

Ae = (points(1,2) - points(1,1))*(points(2,3) - points(2,2)); 

zeta = [-1/sqrt(3) 1/sqrt(3)];
eta = [-1/sqrt(3) 1/sqrt(3)];

w = [1 1];

x_val = 3*[0.2113 0.7887];

y_val = 2*[0.2113 0.7887]-1;

J_det = 1.5;

K = zeros(8);

for p = 1:2
    for q = 1:2

        x = x_val(p);
        y = y_val(q);
        
        H = (1/Ae)*[(y - y_e(4)), 0, -(y - y_e(4)), 0 (y - y_e(1)), 0 , -(y - y_e(1)), 0;
                    0, (x - x_e(2)), 0, -(x - x_e(1)), 0, (x - x_e(1)), 0, -(x - x_e(2));
                    (x - x_e(2)) (y - y_e(4)) -(x - x_e(1)) -(y - y_e(4)) (x - x_e(1)) (y - y_e(1)) -(x - x_e(2)) -(y - y_e(1))];

        K_temp = w(p)*w(q)*J_det*transpose(H)*D*H;
        
        K = K + K_temp;

    end
end

% Value of Element Stiffness Matrix

K

