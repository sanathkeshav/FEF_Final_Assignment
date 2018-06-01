%%% Problem 3
clear; close all; clc

addpath('Func_ReferenceElement')

mu = 1000;

elemV = 0; degreeV = 2; degreeP = 1;
if elemV == 11
    elemP = 1;
else
    elemP = elemV;
end
referenceElement = SetReferenceElementStokes(elemV,degreeV,elemP,degreeP);

Nr = 20;
Ntheta = 15;
[X,T,XP,TP,dbc1,dbc2] = createMesh(Nr,Ntheta,referenceElement);

%figure; plotMesh(T,X,elemV,'b-');

[Conv] = ConvVel(X);
DF = 5;
sig_F = 0.25;
DG = 15;
sig_G = 2;
sig_GF = 0.5;
dvalue = 80;

ngaus=9;
[pospg,wpg] = Quadrature(elemV,ngaus);
[N,Nxi,Neta] = ShapeFunc(elemV,degreeV,pospg);

M = CreMassMat(X,T,pospg,wpg,N,Nxi,Neta);
K = CreStiffMat(X,T,Conv,pospg,wpg,N,Nxi,Neta);

numnp = size(X,1); 
boundary_vals = 0*ones(size(dbc1,2),1);
internal_dofs = setdiff(1:numnp,dbc1);

disp(' ')
t_end = input('End time (default: 0.5) =  ');
nstep = input('Number of time steps (default: 50)= ');
if isempty(t_end)
    t_end = 0.5;
end
if isempty(nstep)
    nstep = 50;
end

dt = t_end/nstep;

% INITIAL CONDITION: Zero
u = zeros(numnp,nstep+1);
u(dbc1,:) = dvalue;
ug = zeros(numnp,nstep+1);


% Matrices arising from the discretization
[Ku,Gu,fu] = StokesSystem2(X,T,XP,TP,referenceElement);
Ku = mu*Ku;
[ndofP,ndofV] = size(Gu);
Xe_ref = referenceElement.Xe_ref;
[Tf,Tu] = boundaryMatrices(X,T,elemV,degreeV,Xe_ref);

[A_DirBC, b_DirBC, nDir, confined] = BC(X,dbc1,dbc2,ndofV);

velo = zeros(size(X,1),2,nstep);

tol = 1e-10;
for n =1:nstep
    error_v = 10;
    error_f = 10;
    z = 0;
    f_aux = u(:,n);
    velo_aux_prev = zeros(size(X,1),2);
    f_aux_prev = u(:,n);
    while error_v > tol || error_f > tol
        Atot = [-Ku+Tu          A_DirBC'
            A_DirBC    zeros(nDir,nDir)];
        btot = [fu - Tf*f_aux ; b_DirBC];
        sol = Atot\btot;
        velo_aux = reshape(sol(1:ndofV), 2, [])';
        Conv = velo_aux;
        C = CreConvMat(X,T,Conv,pospg,wpg,N,Nxi,Neta);
        aux = zeros(numnp,1);
        A = M + (dt/2)*C + (dt/2)*DF*K + (dt/2)*sig_F*M ;
        B = dt*C - dt*DF*K - dt*sig_F*M;
        F = B*u(:,n);
        F = F - A(:, dbc1) * boundary_vals;
        internal_dofs = setdiff(1:numnp,dbc1);
        aux(internal_dofs) = A(internal_dofs, internal_dofs) \ F(internal_dofs); % Solve
        aux(dbc1) = boundary_vals;
        f_aux = u(:,n) + aux(1:numnp);
        error_v = max(max(abs(velo_aux - velo_aux_prev)));
        error_f = max(abs(f_aux - f_aux_prev));
        velo_aux_prev = velo_aux;
        f_aux_prev = f_aux;
        z = z + 1
    end
    velo(:,:,n+1) = velo_aux;
    u(:,n+1) = f_aux;
    A1 = M + (dt/2)*DG*K + (dt/2)*sig_G*M;
    A2 = (dt/2)*sig_GF*M;
    B1 = -dt*DG*K - dt*sig_G*M;
    B2 = dt*sig_GF*M;
    internal_dofs = 1:numnp;
    F = A2*(u(:,n+1) - u(:,n)) + B1*ug(:,n) + B2*u(:,n);
    aux = A1 \ F;
    ug(:,n+1) = ug(:,n) + aux;
    
end

velort = zeros(size(X,1),2);
for i = 1:size(X,1)
    rot = [X(i,1)/sqrt(X(i,1)^2 + X(i,2)^2) X(i,2)/sqrt(X(i,1)^2 + X(i,2)^2)...
        ; -X(i,2)/sqrt(X(i,1)^2 + X(i,2)^2) X(i,1)/sqrt(X(i,1)^2 + X(i,2)^2)];
    velort(i,:) = (rot*velo(i,:,end)')';
end

%%% Postprocess

%figure(1)
nPt = size(X,1);
quiver(X(1:nPt,1),X(1:nPt,2),velo(1:nPt,1,end),velo(1:nPt,2,end));
axis equal; axis tight

%figure(2)
PlotResults(X,T,velo(:,1,end),referenceElement.elemV,referenceElement.degreeV)
%figure(3)
PlotResults(X,T,velo(:,2,end),referenceElement.elemV,referenceElement.degreeV)
v = sqrt((velo(:,1,end).^2)+(velo(:,2,end).^2));
%figure(4)
PlotResults(X,T,v,referenceElement.elemV,referenceElement.degreeV)

if degreeV == 2
    Nr = 2*Nr;
    Ntheta = 2*Ntheta;
end
%figure(5); clf;
[xx,yy,sol] = MatSol(X,Ntheta,Nr,u(:,end));
surface(xx,yy,sol);
view([40,30])
axis auto
grid on;

%figure(11); clf;
[xx,yy,sol] = MatSol(X,Ntheta,Nr,ug(:,end));
surface(xx,yy,sol);
view([40,30])
grid on;

figure(10); clf;
    title("F evolution"); axis equal tight;
    peli = moviein(nstep+1);
    axis auto
    grid on;
    for n=1:nstep+1
        [xx,yy,sol] = MatSol(X,Ntheta,Nr,u(:,n));
        surf(xx,yy,sol);
        pause(0.1)
        peli(:,n) = getframe;
    end
    
 figure(11); clf;
    title("G evolution"); axis equal tight;
    peli = moviein(nstep+1);
    axis auto
    grid on;
    for n=1:nstep+1
        [xx,yy,sol] = MatSol(X,Ntheta,Nr,ug(:,n));
        surf(xx,yy,sol);
        pause(0.1)
        peli(:,n) = getframe;
    end

PlotResults(X,T,velort(:,1),referenceElement.elemV,referenceElement.degreeV)
PlotResults(X,T,velort(:,2),referenceElement.elemV,referenceElement.degreeV)