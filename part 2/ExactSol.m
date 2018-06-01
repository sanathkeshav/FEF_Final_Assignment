function [u,v,u_x,u_y,v_x,v_y,p] = ExactSol(pt)
x = pt(:,1); 
y = pt(:,2); 


u = 2*x.^2.*(1-x).^2.*y.*(y-1).*(2*y-1); 
u_x = 4*x.*y.*(2*y-1).*(y-1).*(2*x-1).*(x-1); 
u_y = 2*x.^2.*(6*y.^2-6*y+1).*(x-1).^2; 
v = -2*y.^2.*(1-y).^2.*x.*(x-1).*(2*x-1); 
v_x = -2*y.^2.*(y-1).^2.*(6*x.^2-6*x+1); 
v_y = -4*x.*y.*(2*y-1).*(y-1).*(2*x-1).*(x-1); 
p = x.*(1-x);