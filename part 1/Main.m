clear;
close all;
clc;

disp(' ')
Nr = input('Number of elements on R direction (default: 20)= ');
Ntheta = input('Number of elements on theta direction (default: 10) = ');
if isempty(Nr)
    Nr = 20;
end
if isempty(Ntheta)
    Ntheta = 10;
end
[X,T,dbc] = createMesh(Nr,Ntheta);
numnp = size(X,1); 

% CONVECTION VELOCITY
[Conv] = ConvVel(X); 
% Parameters
DF = 5;
sig_F = 0.25;
DG = 15;
sig_G = 2;
sig_GF = 0.5;
dvalue = 80;


% TIME-INTEGRATION SCHEME
% Number of Gauss points (numerical quadrature)
ngaus=4;
% Quadrature
[pospg,wpg] = Quadrature(ngaus);
% Shape Functions
[N,Nxi,Neta] = ShapeFunc(pospg);

% COMPUTATION OF THE MATRICES 
%disp ('Computation of mass matrix M, obtained by discretizing the term (w,u)')
M = CreMassMat(X,T,pospg,wpg,N,Nxi,Neta);
%disp ('Computation of convection matrix C, obtained by discretizing the term (a�grad(w),u)')
C = CreConvMat(X,T,Conv,pospg,wpg,N,Nxi,Neta);
%disp ('Computation of stiffness matrix K, obtained by discretizing the term (a�grad(w),a�grad(u))')
K = CreStiffMat(X,T,Conv,pospg,wpg,N,Nxi,Neta);

% Vectors related to source term
v1 = CreVect1(X,T,pospg,wpg,N,Nxi,Neta);          % (w,s) 
v2 = CreVect2(X,T,Conv,pospg,wpg,N,Nxi,Neta);     % (a�grad(w),s) 


% DATA FOR THE TRANSIENT ANALYSIS
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
u(dbc,:) = dvalue;
ug = zeros(numnp,nstep+1);

% COMPUTATION OF MATRICES FOR THE TRANSIENT ANALYSIS
 
A = M + (dt/2)*C + (dt/2)*DF*K + (dt/2)*sig_F*M ;
B = dt*C - dt*DF*K - dt*sig_F*M;

A1 = M + (dt/2)*DG*K + (dt/2)*sig_G*M;
A2 = (dt/2)*sig_GF*M;
B1 = -dt*DG*K - dt*sig_G*M;
B2 = dt*sig_GF*M;

for n = 1:nstep
    aux = zeros(numnp,1);
    F = B*u(:,n);
    boundary_vals = 0*ones(size(dbc,2),1);
    internal_dofs = setdiff(1:numnp,dbc);
    F = F - A(:, dbc) * boundary_vals;
    aux(internal_dofs) = A(internal_dofs, internal_dofs) \ F(internal_dofs); % Solve
    aux(dbc) = boundary_vals;
    u(:,n+1) = u(:,n) + aux(1:numnp);
    
    aux = zeros(numnp,1);
    F = A2*(u(:,n+1)-u(:,n)) + B1*ug(:,n) + B2*u(:,n);
    internal_dofs = 1:numnp;
    aux(internal_dofs) = A1(internal_dofs, internal_dofs) \ F(internal_dofs);
    ug(:,n+1) = ug(:,n) + aux(1:numnp);
end


% POSTPROCESS
if max(u(:,nstep+1))<100 &  min(u(:,nstep+1))>-100
        % Solution at time t=t_end
    	figure(1); clf;
    	[xx,yy,sol] = MatSol(X,Ntheta,Nr,u(:,nstep+1));
        surface(xx,yy,sol);
    	view([40,30])
    	grid on;
        
        figure(3); clf;
    	[xx,yy,sol] = MatSol(X,Ntheta,Nr,ug(:,nstep+1));
        surface(xx,yy,sol);
    	view([40,30])
    	grid on;
    
        % Contour plot of the solution at time t = t_end
    	figure(2); clf;
        set(gca,'FontSize',12);
        [xx,yy,sol] = MatSol(X,Ntheta,Nr,u(:,nstep+1));
    	[C,h]=contour(xx,yy,sol,15);
    	clabel(C,h);
     	axis equal tight;
        
        figure(4); clf;
        set(gca,'FontSize',12);
        [xx,yy,sol] = MatSol(X,Ntheta,Nr,ug(:,nstep+1));
    	[C,h]=contour(xx,yy,sol,15);
    	clabel(C,h);
    	axis equal tight;
    
    % Movie
    figure(10); clf;
    title("F evolution");  axis equal tight;
    peli = moviein(nstep+1);
    axis auto
    grid on;
  
    for n=1:nstep+1
        [xx,yy,sol] = MatSol(X,Ntheta,Nr,u(:,n));
        surf(xx,yy,sol);
        pause(0.1)
        peli(:,n) = getframe;
    end
    % Movie
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
end

nPt = size(X,1);
quiver(X(1:nPt,1),X(1:nPt,2),Conv(1:nPt,1),Conv(1:nPt,2));
axis equal; axis tight


