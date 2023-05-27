% Assignment 6 Q1
clear
load("ps_06_data.mat");

xvals = (1:50);
yvals = -1*(1:50);

[X,Y] = meshgrid(xvals,yvals);

contourf(X,Y,rho,100,'edgecolor','none');
colorbar();