% Assignment 6 Q2
clear
load("ps_06_data.mat")
gamma = 0.2;
E = 70e+3; %[in MPa]
nu = 0.3;
volfrac = 0.3;
penal = 3;

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


adj_stress_vm = stress_vm./(rho.^penal);
rel_stress_vm = zeros(ely,elx);
for elx =1:nelx
    for ely = 1:nely
        rel_stress_vm(ely,elx) = stress_vm(ely,elx)/((rho(ely,elx)^penal)*(1-gamma+(gamma/rho(ely,elx))));
    end
end

%contourf(X,Y,stress_vm,100,'edgecolor','none')
figure(1)
contourf(X,Y,adj_stress_vm,100,'edgecolor','none')
title('Adjusted Von-Misses Stress')
colorbar()

figure(2)
contourf(X,Y,rel_stress_vm,100,'edgecolor','none')
title('Relaxed Von-Mises Stress')
colorbar()