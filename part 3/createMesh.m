function [X,T,XP,TP,dbc1,dbc2] = createMesh(Nr,Ntheta,referenceElement)
%function [XR,TR,dbc] = createMesh(Nr,Ntheta,referenceElement)
%Nr = 2 ;
%Ntheta = 2 ;

elemV = referenceElement.elemV;
degreeV = referenceElement.degreeV;
elemP = referenceElement.elemP;
degreeP = referenceElement.degreeP;

if degreeP == 1
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
            ele_data(ind_ele,:) = [ind_ele , Node_number(ii,jj) , Node_number(ii+1,jj), Node_number(ii+1,jj+1), Node_number(ii,jj+1)] ;
            %ele_data(ind_ele,:) = [ind_ele , Node_number(ii+1,jj) , Node_number(ii+1,jj+1), Node_number(ii,jj+1), Node_number(ii,jj)] ;
        end
    end
    
    % filename = 'elements_trial.txt' ;
    % dlmwrite(filename,ele_data,'\t') ;
    % eval(['open ' filename])
    
    XR = node_data(:,2:3) ;
    TR = ele_data(:,2:end) ;
elseif degreeP == 2
    Nr = Nr*2;
    Ntheta = Ntheta*2;
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
    for jj = 1:2:Nr-1
        for ii = 1:2:Ntheta-1
            ind_ele = ind_ele +1 ;
            ele_data(ind_ele,:) = [ind_ele , Node_number(ii+2,jj) , Node_number(ii+2,jj+2),...
                Node_number(ii,jj+2), Node_number(ii,jj), Node_number(ii+2,jj+1), Node_number(ii+1,jj+2),...
                Node_number(ii,jj+1), Node_number(ii+1,jj), Node_number(ii+1,jj+1)] ;
        end
    end
    
    % filename = 'elements_trial.txt' ;
    % dlmwrite(filename,ele_data,'\t') ;
    % eval(['open ' filename])
    
    XR = node_data(:,2:3) ;
    TR = ele_data(:,2:end) ;
end
XP = XR;
TP = TR;

clear ele_data;
clear node_data;
if degreeV == 1
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
elseif degreeV == 2
    Nr = Nr*2;
    Ntheta = Ntheta*2;
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
    for jj = 1:2:Nr-1
        for ii = 1:2:Ntheta-1
            ind_ele = ind_ele +1 ;
            ele_data(ind_ele,:) = [ind_ele , Node_number(ii+2,jj) , Node_number(ii+2,jj+2),...
                Node_number(ii,jj+2), Node_number(ii,jj), Node_number(ii+2,jj+1), Node_number(ii+1,jj+2),...
                Node_number(ii,jj+1), Node_number(ii+1,jj), Node_number(ii+1,jj+1)] ;
        end
    end
    
    % filename = 'elements_trial.txt' ;
    % dlmwrite(filename,ele_data,'\t') ;
    % eval(['open ' filename])
    
    XR = node_data(:,2:3) ;
    TR = ele_data(:,2:end) ;
end
X = XR;
T = TR;
%plotMesh(TR,XR,0,'r',1)

dbc1 = [];
for j = 1:size(XR,1)
    rxy = sqrt(XR(j,1)^2 + XR(j,2)^2);
    if rxy <= 25.001 && rxy>= 24.999
        dbc1 = [dbc1 j];
    end
end

dbc2 = [];
for j = 1:size(XR,1)
    rxy = sqrt(XR(j,1)^2 + XR(j,2)^2);
    if rxy <= 15.001 && rxy>= 14.999
        dbc2 = [dbc2 j];
    end
end

end