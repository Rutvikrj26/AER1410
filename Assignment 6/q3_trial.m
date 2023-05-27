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
total_strain_energy = 0;
stress_vm = zeros(nelx,nely);
Ae = 1;
dc = zeros(nely,nelx);
for elx = 1:nelx
    for ely = 1:nely
        disp = [ dx(ely+1,elx) dy(ely+1,elx); dx(ely+1,elx+1) dy(ely+1,elx+1); dx(ely,elx+1) dy(ely,elx+1); dx(ely,elx) dy(ely,elx)];
        D = (rho(ely,elx)^penal)*(E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2];
        strain = 0.5*[disp(2,1)+disp(3,1)-disp(4,1)-disp(1,1);
                      disp(3,2)+disp(4,2)-disp(1,2)-disp(2,2);
                      disp(2,2)+disp(3,1)+disp(3,2)+disp(4,1)-disp(1,1)-disp(1,2)-disp(2,1)-disp(4,2)];
        stress_temp = D*strain;
        strain_energy = 0.5*stress_temp'*strain;
        total_strain_energy = total_strain_energy + volfrac*strain_energy;
        dc(ely,elx) = -penal*rho(ely,elx)^(penal-1)*2*strain_energy;
    end
end

hold on
limit = [-0.3,-2];
contourf(X,Y,dc,100,'edgecolor','none')
colorbar()
contour(X,Y,dc,limit,'edgecolor','black')
hold off


%% Comment on the Distribution of Gradient
% As can be clearly seen in the Contour plot that the gradient has an area
% that has less negative values than others. This means, that part can be made more
% stiff by adding very less amount of material, and thus more material will
% be added there.
% The contour has an area enclosed by black thick lines, inside which more
% material will be added, and at the other locations, material will be
% removed.