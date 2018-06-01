function [A, b2, nDir, confined] = BC(X,dbc1,dbc2,n)

confined = 1; 
nodesDy0 = dbc2';
nodesDy1 = dbc1';
nodesDirBC = [nodesDy0; nodesDy1];

nDir = 2*length(nodesDirBC);
C = [2*nodesDirBC - 1; 2*nodesDirBC];
C1=reshape(C,nDir/2,2);
C2 =reshape(C1',nDir,1);
A = zeros(nDir,n); 
A(:,C2) = eye(nDir); 

vx0 = (-0.15)* X(dbc2(:),1)./sqrt(X(dbc2(:),1).^2 + X(dbc2(:),2).^2);
vy0 = (-0.15)* X(dbc2(:),2)./sqrt(X(dbc2(:),1).^2 + X(dbc2(:),2).^2);
vx1 = (-0.3)* X(dbc1(:),1)./sqrt(X(dbc1(:),1).^2 + X(dbc1(:),2).^2);
vy1 = (-0.3)* X(dbc1(:),2)./sqrt(X(dbc1(:),1).^2 + X(dbc1(:),2).^2);

% for i = 1:size(dbc1,2)
%     if X(dbc1(i),1) <= 0
%         vx1(i) = -(-0.3)* X(dbc1(i),1)/sqrt(X(dbc1(i),1)^2 + X(dbc1(i),2)^2);
%         vy1(i) = -(-0.3)* X(dbc1(i),2)/sqrt(X(dbc1(i),1)^2 + X(dbc1(i),2)^2);
%     else
%         vx1(i) = (-0.3)* X(dbc1(i),1)/sqrt(X(dbc1(i),1)^2 + X(dbc1(i),2)^2);
%         vy1(i) = (-0.3)* X(dbc1(i),2)/sqrt(X(dbc1(i),1)^2 + X(dbc1(i),2)^2);
%     end
% end
% 
% for i = 1:size(dbc2,2)
%     if X(dbc2(i),1) <= 0
%         vx0(i) = -(-0.15)* X(dbc2(i),1)/sqrt(X(dbc2(i),1)^2 + X(dbc2(i),2)^2);
%         vy0(i) = -(-0.15)* X(dbc2(i),2)/sqrt(X(dbc2(i),1)^2 + X(dbc2(i),2)^2);
%     else
%         vx0(i) = (-0.15)* X(dbc2(i),1)/sqrt(X(dbc2(i),1)^2 + X(dbc2(i),2)^2);
%         vy0(i) = (-0.15)* X(dbc2(i),2)/sqrt(X(dbc2(i),1)^2 + X(dbc2(i),2)^2);
%     end
% end

b1 = [vx0 vy0;...
    vx1   vy1];
b2 =reshape(b1',nDir,1);


