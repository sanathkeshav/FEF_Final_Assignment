function [Conv] = ConvVel(X)
% Conv = ConvVel(X)
% Given a matrix of nodal coordinates X,
% this functions provides a nodal velocity field Conv
%
Conv = zeros(size(X));
Conv(:,1) = -(1/1500)*X(:,1).*sqrt(X(:,1).^2+X(:,2).^2);
Conv(:,2) = -(1/1500)*X(:,2).*sqrt(X(:,1).^2+X(:,2).^2);
end

