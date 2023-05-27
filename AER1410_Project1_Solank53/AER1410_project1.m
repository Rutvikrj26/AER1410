%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
clear
%% Basic Inputs

% Student Number = 1008247531;

A = 5; B = 3; C = 1;

%% Setting up the Geometry

base_width = 50 + 70 + 10*B;
tot_height = 250 + 10*A;

P1 = 3;
P2 = P1/(0.5*(C+1));

nelx = base_width;
nely = tot_height;
volfrac = 0.133;
penal = 3;
rmin = 6;
ft = 1;

%% MATERIAL PROPERTIES
E0 = 69;
Emin = 1e-9;
nu = 0.3;

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
F(2*(nely+1)) = -P1;
F(2*(nely+1)*(nelx+1)) = -P2;
fixeddofs   = union(2*[(nely+1)*(0)+1]:2*(nely+1):2*[(nely+1)*(nelx)+1], 2*[(nely+1)*(0)+1]-1:2*(nely+1):2*[(nely+1)*(nelx)+1]-1);
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
for ely = 1:nely
    for elx = 1:nelx
        if (ely < tot_height-50) && (elx < 50 || elx > 100)
            passive(ely,elx) = 1;
        else
            passive(ely,elx) = 0;
        end
    end
end

x(find(passive))=0.001;

xPhys = x;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    xnew(find(passive)) = 0.001;
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES 
  fig = colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  print(['Frame ' num2str(loop)],'-dpng');
end


[Y,I] = sort(x(:), 'descend');
vt = floor(((volfrac - 0.001) * nelx * nely) / (1 - 0.001));
xd(I(1:vt)) = 1;
xd(I(vt+1:end)) = 0.001;
xd = reshape(xd, nely, nelx);

% sK = reshape(KE(:)*(Emin+xd(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
% K = sparse(iK,jK,sK); K = (K+K')/2;
% U(freedofs) = K(freedofs,freedofs)\F(freedofs);

colormap(gray); imagesc(1-xd); caxis([0 1]); axis equal; axis off; drawnow;
print(['Frame ' num2str(loop+1)],'-dpng');

%% Now time to generate the animated gif
GifName = 'project_1.gif';
delay = 0.5;    % Delay between frames (s)
for ii = 1:loop+1
    [A, ~] = imread(['Frame ' num2str(ii) '.png']);
    [X, map] = rgb2ind(A, 256);
    if ii == 1
        imwrite(X, map, GifName, 'gif', 'LoopCount', inf, 'DelayTime', delay)
    else
        imwrite(X, map, GifName, 'gif', 'WriteMode', 'append', 'DelayTime', delay)
    end
end

