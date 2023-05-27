% Assignment 5 Q2
clear
l1 = 0;l2 = 1e+5;
lmid = 0.5*(l1+l2);

f = @(x) exp(x/25) * (sin(x*pi/48) +2);

while abs(l1-l2)>1e-4
lmid = 0.5*(l1+l2);
%f(lmid)
    if f(lmid) > 25
        l2 = lmid;
    else
        l1 = lmid;
    end
end

disp('Final Root of the Equation:')

lmid

