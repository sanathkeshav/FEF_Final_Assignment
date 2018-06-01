% This program solves Stokes problem in a square domain


clear; close all; clc

addpath('Func_ReferenceElement')

dom = [0,1,0,1]; 
mu = 1000; 

% Element type and interpolation degree
% (0: quadrilaterals, 1: triangles, 11: triangles with bubble function)
elemV = 0; degreeV = 2; degreeP = 1;
% elemV = 1; degreeV = 2; degreeP = 1;
% elemV = 11; degreeV = 1;  degreeP = 1; 
if elemV == 11
    elemP = 1; 
else
    elemP = elemV; 
end
referenceElement = SetReferenceElementStokes(elemV,degreeV,elemP,degreeP); 

%Nr = cinput('Number of elements in each direction',10);
Nr = 30;
Ntheta = 20; 
[X,T,XP,TP,dbc1,dbc2] = createMesh(Nr,Ntheta,referenceElement);

figure; plotMesh(T,X,elemV,'b-');
figure; plotMesh(TP,XP,elemP,'r-');

% Matrices arising from the discretization
[K,G,f] = StokesSystem2(X,T,XP,TP,referenceElement);
% [K,G,f] = StokesSystem2(X,T,XP,TP,referenceElement);
K = mu*K; 
[ndofP,ndofV] = size(G); 

% Matrix and r.h.s vector to impose Dirichlet boundary conditions using
% Lagrange multipliers
[A_DirBC, b_DirBC, nDir, confined] = BC(X,dbc1,dbc2,ndofV);

% Total system of equations
if confined
   nunkP = ndofP-1;
   disp(' ')
   disp('Confined flow. Pressure on lower left corner is set to zero');
   G(1,:) = [];
else
   nunkP = ndofP;
end
Atot = [K          A_DirBC'             G'
        A_DirBC    zeros(nDir,nDir)     zeros(nDir,nunkP)
        G          zeros(nunkP,nDir)    zeros(nunkP,nunkP)];
btot = [f ; b_DirBC ; zeros(nunkP,1)];

sol = Atot\btot; 

% Postprocess
velo = reshape(sol(1:ndofV), 2, [])'; 
if confined
    pres = [0; sol(ndofV+nDir+1:ndofV+nDir+nunkP)]; 
else
	pres = sol(ndofV+nDir+1:ndofV+nDir+nunkP); 
end

velort = zeros(size(X,1),2);
for i = 1:size(X,1)
    %teta = atan(X(i,2)/X(i,1));
    rot = [X(i,1)/sqrt(X(i,1)^2 + X(i,2)^2) X(i,2)/sqrt(X(i,1)^2 + X(i,2)^2)...
        ; -X(i,2)/sqrt(X(i,1)^2 + X(i,2)^2) X(i,1)/sqrt(X(i,1)^2 + X(i,2)^2)];
    velort(i,:) = (rot*velo(i,:)')';
end

% nPt = size(X,1); 
% figure; 
% quiver(X(1:nPt,1),X(1:nPt,2),velo(1:nPt,1),velo(1:nPt,2));
% hold on 
% plot(dom([1,2,2,1,1]),dom([3,3,4,4,3]),'k')
% axis equal; axis tight
% 
% PlotStreamlines(X,velo,dom); 
% 
% PlotResults(X,T,velo(:,1),referenceElement.elemV,referenceElement.degreeV)
% 
% PlotResults(X,T,velo(:,2),referenceElement.elemV,referenceElement.degreeV)
% 
% if degreeP == 0
%     PlotResults(X,T,pres,referenceElement.elemP,referenceElement.degreeP)
% else
%     PlotResults(XP,TP,pres,referenceElement.elemP,referenceElement.degreeP)
% end

nPt = size(X,1); 
figure; 
quiver(X(1:nPt,1),X(1:nPt,2),velo(1:nPt,1),velo(1:nPt,2));
hold on 
% plot(dom([1,2,2,1,1]),dom([3,3,4,4,3]),'k')
axis equal; axis tight

%PlotStreamlines(X,velo,dom); 

PlotResults(X,T,velo(:,1),referenceElement.elemV,referenceElement.degreeV)

PlotResults(X,T,velo(:,2),referenceElement.elemV,referenceElement.degreeV)

%PlotResults(X,T,velort(:,1),referenceElement.elemV,referenceElement.degreeV)
%PlotResults(X,T,velort(:,2),referenceElement.elemV,referenceElement.degreeV)

if degreeP == 0
    PlotResults(XP,TP,pres,referenceElement.elemP,referenceElement.degreeP)
else
    PlotResults(XP,TP,pres,referenceElement.elemP,referenceElement.degreeP)
end
