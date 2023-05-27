% Assignment-4 Q-2

% [ all units in GPa, mm, KN]

clear
  nelx = 2;
  nely = 2;

  x = [0.45 0.7;0.35 0.6];
  penal = 3;
  ele_X_coord = [0 10 10 0; 10 20 20 10; 0 10 10 0; 10 20 20 10];
  ele_Y_coord = [0 0 10 10; 0 0 10 10; 10 10 20 20; 10 10 20 20];

  ele_nodes = [1 2 5 4; 2 3 6 5; 4 5 8 7; 5 6 9 8];

  [U]=[0.01 0 0.015 0.003 0.02 0.008 0.005 0.004 0.015 0.006 0.03 0.009 0.015 0.004 0.025 0.005 0.035 0.006];         
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  c = 0.;
  curr_ele = 1;
  for ely = 1:nely
    for elx = 1:nelx
      curr_nodes = ele_nodes(curr_ele,:);
      Ue = [];
      for node = curr_nodes
        Ue = [Ue; U(2*node-1); U(2*node)];
      end
      KE = stiffness(ele_X_coord(curr_ele,:), ele_Y_coord(curr_ele,:));
      c = c + x(ely,elx)^penal*Ue'*KE*Ue;
      dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
      curr_ele = curr_ele + 1;
    end
  end

dc

function K = stiffness(X,Y)
E = 72; % [in GPa]
nu = 0.3;
D = (E/((1+nu)*(1-2*nu)))*[1-nu nu 0; nu 1-nu 0;0 0 0.5-nu];
coord = [X',Y'];
K = zeros(8,8);
for i = 1:2
for j = 1:2
eta = (2*i-3)/sqrt(3);
zeta = (2*j-3)/sqrt(3);
J = (1/4)*[eta-1 1-eta 1+eta -eta-1; zeta-1 -zeta-1 1+zeta 1-zeta]*coord;
H = (1/4)*[eta-1 1-eta 1+eta -eta-1; zeta-1 -zeta-1 1+zeta 1-zeta];
H = J\H;
H = [H(1,1) 0 H(1,2) 0 H(1,3) 0 H(1,4) 0; 0 H(2,1) 0 H(2,2) 0 H(2,3) 0 H(2,4); H(2,1) H(1,1) H(2,2) H(1,2) H(2,3) H(1,3) H(2,4) H(1,4)];
K = K + det(J)*H'*D*H;
end
end
end
 
