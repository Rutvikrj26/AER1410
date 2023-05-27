% Plotting the Solution
n = linspace(-10,20,100);
r = linspace(-10,20,100);
[X,Y] = meshgrid(n,r);
Z = X.^2 - 8.*X + Y.^2 - 4.*Y + 15;

contour(X,Y,Z)
colorbar()
xlabel('n')
ylabel('r')