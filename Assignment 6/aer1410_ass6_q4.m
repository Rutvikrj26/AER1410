% Assignment 6 Q4

% Assignment 6 Q2
clear
load("ps_06_data.mat")
gamma = 0.2;
E = 70e+3;
nu = 0.3;
volfrac = 0.3;
penal = 3;
sigma_y = 325;

[nelx, nely] = size(rho);

xvals = [1:nelx] - 0.5;
yvals = -xvals;

[X,Y] = meshgrid(xvals,yvals);

stress_vm = zeros(nelx,nely);
Ae = 1;
for elx = 1:nelx
    for ely = 1:nely
        disp = [ dx(ely+1,elx) dy(ely+1,elx); dx(ely+1,elx+1) dy(ely+1,elx+1); dx(ely,elx+1) dy(ely,elx+1); dx(ely,elx) dy(ely,elx)];
        D = (rho(ely,elx)^penal)*(E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2];
        strain = 0.5*[disp(2,1)+disp(3,1)-disp(4,1)-disp(1,1);
                      disp(3,2)+disp(4,2)-disp(1,2)-disp(2,2);
                      disp(2,2)+disp(3,1)+disp(3,2)+disp(4,1)-disp(1,1)-disp(1,2)-disp(2,1)-disp(4,2)];
        stress_temp = D*strain;
        stress_xx = stress_temp(1);
        stress_yy = stress_temp(2);
        stress_xy = stress_temp(3);
        stress_vm(ely,elx) = sqrt((stress_xx.^2 + stress_yy.^2) + 3*stress_xy.^2 - stress_xx.*stress_yy);        
    end
end
P_norm = [];
for P = 3:500
P_norm = [P_norm (sum(sum((max(0,stress_vm)).^P)))^(1/P)];
end

plot(3:500,P_norm)

%% Comments on using Von Mises Stress for Stress Limits
% Using Von-Misses stress directly will lead to a number of problems.
% Von-misses stress is the equivallent stress calculated after resolving
% the stresses in all three directions. The idea behing using the P-Norms
% is to calculate the one global bounding limit while applying the stress
% limits. In case we use the Von-Mises Stess for calculating the P-Norm and
% P-Mean values, it will lead to a problem of convergence.
%
% It can be clearly seen that plot does not show up after 130 eventhough it is supposed to show up till 500(according to code), because the
% values after 130 become infinity. This is the key computational issue in
% using the Von-Mises stress directly.