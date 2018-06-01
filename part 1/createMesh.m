function [XR,TR,dbc] = createMesh(Nr,Ntheta) 

 %Nr = 20 ;
 %Ntheta = 10 ;

r = linspace(15,25,Nr+1) ;
theta_max = 10*pi/180 ;
theta = linspace(-theta_max,theta_max,Ntheta+1) ;

[R,TH ] = meshgrid(r,theta) ;
X = R.*sin(TH) ;
Y = R.*cos(TH) ;

Node_number = 1:numel(X) ;
Node_number = reshape(Node_number,size(X)) ;

node_data = zeros(numel(X),4) ;
filename = 'nodes_trial.txt' ;
for i_node = 1:numel(X)
	node_data(i_node,:) = [Node_number(i_node) , X(i_node) , Y(i_node), 0] ;
end

% dlmwrite(filename,node_data,'\t')
% eval(['open ' filename])

ind_ele = 0 ; 
for jj = 1:Nr
	for ii = 1:Ntheta
		ind_ele = ind_ele +1 ;
		ele_data(ind_ele,:) = [ind_ele , Node_number(ii+1,jj) , Node_number(ii+1,jj+1), Node_number(ii,jj+1), Node_number(ii,jj)] ;
	end
end

% filename = 'elements_trial.txt' ;
% dlmwrite(filename,ele_data,'\t') ;
% eval(['open ' filename])

XR = node_data(:,2:3) ;
TR = ele_data(:,2:end) ;
plotMesh(TR,XR,0,'b',0)

dbc = [];
for j = 1:size(XR,1)
    rxy = sqrt(XR(j,1)^2 + XR(j,2)^2);
    if rxy <= 25.001 && rxy>= 24.999
        dbc = [dbc j];
    end
end

end