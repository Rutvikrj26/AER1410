% AER1410 Assignment-8 Q-2


syms x
fun = sin(pi*x/2) + 20/(x+2) + exp(-x)*(x^4) + (x/10)^5;
p = 3;
error = 2;
while error >= 0.1
p = p +1;
interval = 10/(p-1);
points = linspace(0,10+interval,p+1);
psi = 0;
for i = 1:p
    yp = subs(fun,points(i));
    psi = psi + max(yp - (yp/interval)*abs(x-points(i)),0);
end

error = vpaintegral((fun-psi)^2,[0 10]);
end

disp("Number of Basis Function Needed :")
disp(p)

disp("error value :")
disp(error)

disp("Values of the Basis Function :")
disp(subs(fun,points'))

vals = linspace(0,10,100);
plot(vals,subs(fun,vals))
hold on;
plot(vals,subs(psi,x,vals))
hold off