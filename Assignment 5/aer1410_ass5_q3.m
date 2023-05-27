% Assignment 5 Q3
clear

% The provided data is loaded into the script via this command. Please make
% sure it is present with right name in the directory while running this script. 

load('problem_set_05_data.mat') 

volfrac = 0.3; nel = 20;

l1 = 0; l2 = 100000; move = 0.25;
while (l2-l1 > 1e-4)
lmid = 0.5*(l2+l1);
density_new = max(0.001,max(density-move,min(1.,min(density+move,density.*(-sensitivity./lmid).^0.75))));
if sum(sum(density_new)) - volfrac*nel > 0
l1 = lmid;
else
l2 = lmid;
end
end

final_matrix  = zeros(20,4);

final_matrix(:,1) = density;
final_matrix(:,2) = sensitivity;
final_matrix(:,3) = density_new;
final_matrix(:,4) = density_new - density;

disp('The Final Matrix:')

final_matrix