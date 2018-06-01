function [Tf,Tu] = boundaryMatrices(X,T,elem,degree,Xe_ref)
% [Tf,Tu] = boundaryMatrices(X,T,elem,p,nen,Xe_ref)
% Matrices arising from the discretization of the terms on the boundary. 
% - Discretizing the term div(sigma_m) yields Tf*F, where F is the vector
% of nodal densities [F1,F2,... Fn]  
% - Discretizing the term Tm tields Tu*u, where u is the vector of nodal
% velocities [u1,v1,u2,v2,... un,vn]. Here, u stands for the horizontal
% velocity and v is the vertical one. 
% 
% Warning: note that straight-sided elements are considered
%
% X: matrix of nodal coordinates
% T: connectivity matrix
% elem: type of element (0: quadrilatera, 1: triangles)
% degree: interpolation degree
% Xe_ref: nodal coordinates of the reference element


nen = size(Xe_ref,1); 

% 1D Gauss quadrature
ngaus1D = 3; 
[zgp1D, wgp1D] = Quadrature_1D(ngaus1D); 

nDof = size(X,1); 
nElem = size(T,1); 

nedofV = 2*nen; 
nedofF = nen;
ndofV = 2*nDof; 
ndofF = nDof; 

r = sqrt(X(:,1).^2+X(:,2).^2);
nodesB = find(abs(r-25) < 1e-6);
elemsB = zeros(nElem,3); indB = 1; 
for ielem = 1:nElem
    Te = T(ielem,:);
    Xe = X(Te,:); 
    aux = intersect(Te,nodesB); 
    if length(aux) == 2
        node1 = find(Te == aux(1)); 
        node2 = find(Te == aux(2)); 
        if Xe(node1,1) > Xe(node2,1)
            elemsB(indB,:) = [ielem,node1,node2]; 
        else
            elemsB(indB,:) = [ielem,node2,node1]; 
        end
        indB = indB+1; 
    end
end
elemsB = elemsB(1:indB-1,:); 

n = length(elemsB*nedofV^2+2); 
coef_Tu  = zeros(1,n); indTu_i   = zeros(1,n); indTu_j   = zeros(1,n); 
indTu_i(1) = 1;     indTu_j(1) = 1;     coef_Tu(1) = 0; 
indTu_i(2) = ndofV; indTu_j(2) = ndofV; coef_Tu(2) = 0; 
indTu = 3; 
coef_Tf  = zeros(1,n); indTf_i   = zeros(1,n); indTf_j   = zeros(1,n); 
indTf_i(1) = 1;     indTf_j(1) = 1;     coef_Tf(1) = 0; 
indTf_i(2) = ndofV; indTf_j(2) = ndofF; coef_Tf(2) = 0; 
indTf = 3; 
%Tu = zeros(ndofV,ndofV);
%Tf = zeros(ndofV,ndofF); 

for i = 1:size(elemsB,1) 
    ielem = elemsB(i,1); 
    node1 = elemsB(i,2); 
    node2 = elemsB(i,3); 
    Te = T(ielem,:); 
    Xe = X(Te,:); 
    
    TeV_dof = reshape([2*Te-1; 2*Te],1,nedofV);
    TeF_dof = Te; 
    
    [n,wgp,N,Nxi,Neta] = evaluateAtEdge(node1,node2,Xe,Xe_ref,elem,degree,zgp1D,wgp1D); 
    
    [Tf_e,Tu_e] = EleMat(Xe,n,nedofV,nedofF,ngaus1D,wgp,N,Nxi,Neta);
    
    for irow = 1:nedofV
        for icol = 1:nedofV
            indTu_i(indTu) = TeV_dof(irow);
            indTu_j(indTu) = TeV_dof(icol);
            coef_Tu(indTu) = Tu_e(irow,icol);
            indTu = indTu+1;
        end
        for icol = 1:nedofF
            indTf_i(indTf) = TeV_dof(irow);
            indTf_j(indTf) = TeF_dof(icol);
            coef_Tf(indTf) = Tf_e(irow,icol);
            indTf = indTf+1;
        end
    end
    %Tf(TeV_dof,TeF_dof) = Tf(TeV_dof,TeF_dof) + Tf_e;
    %Tu(TeV_dof,TeV_dof) = Tu(TeV_dof,TeV_dof) + Tu_e;
end
% Create sparse matrices
indTu_i  = indTu_i(1:indTu-1);
indTu_j  = indTu_j(1:indTu-1);
coef_Tu = coef_Tu(1:indTu-1);
Tu = sparse(indTu_i,indTu_j,coef_Tu);
indTf_i  = indTf_i(1:indTf-1);
indTf_j  = indTf_j(1:indTf-1);
coef_Tf = coef_Tf(1:indTf-1);
Tf = sparse(indTf_i,indTf_j,coef_Tf);

Tf = -560*Tf; 
Tu = -0.5*Tu; 



function [n,wgp,N,Nxi,Neta] = evaluateAtEdge(node1,node2,Xe,Xe_ref,elem,degree,zgp1D,wgp1D)

% Normal vector
P1 = Xe(node1,:);  
P2 = Xe(node2,:);
t = P2 - P1; 
h = norm(t);  
n = [t(2), -t(1)]; 
n = n/h; 

% Modified integration weigths 
wgp = wgp1D*h/2; 

% Integration points on the 2D element
P1_ref = Xe_ref(node1,:);  
P2_ref = Xe_ref(node2,:);
zgp = [ 
    (P2_ref(1) - P1_ref(1))/2*zgp1D + (P2_ref(1) + P1_ref(1))/2,...
    (P2_ref(2) - P1_ref(2))/2*zgp1D + (P2_ref(2) + P1_ref(2))/2
    ]; 

% 2D shape functions at the integration points
[N,Nxi,Neta] = ShapeFunc(elem,degree,zgp); 






function [Tf_e,Tu_e] = EleMat(Xe,n,nedofV,nedofF,ngaus,wgp,N,Nxi,Neta)

Tf_e = zeros(nedofV,nedofF);
Tu_e = zeros(nedofV,nedofV);

for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
    Jacob = [
        Nxi_ig*(Xe(:,1))	Nxi_ig*(Xe(:,2))
        Neta_ig*(Xe(:,1))	Neta_ig*(Xe(:,2))
        ];
    res = Jacob\[Nxi_ig;Neta_ig];
    nx = res(1,:);
    ny = res(2,:);
    
    Nn1 = n(1)^2*nx + n(1)*n(2)*ny; 
    Nn2 = n(1)*n(2)*nx + n(2)^2*ny; 
    Nn = reshape([Nn1;Nn2],1,nedofV); 
    Tf_e = Tf_e + Nn'*N_ig*wgp(ig); 
    
	Naux = [reshape([1;0]*N_ig,1,nedofV); reshape([0;1]*N_ig,1,nedofV)];
    Tu_e = Tu_e + Naux'*Naux*wgp(ig); 
end
 

function [z, w] = Quadrature_1D(ngaus)
%
% [z, w] = Quadrature_1D(ngaus)

if ngaus == 1
    z = 0;
    w = 2;
elseif ngaus == 2
    pos1 = 1/sqrt(3);
    z = [-pos1; pos1];
    w = [ 1 1 ];
elseif ngaus == 3
    pos1 = sqrt(3/5);
    z = [-pos1; 0; pos1];
    w = [ 5/9   8/9   5/9 ];
elseif ngaus == 4
    pos1 = sqrt(525+70*sqrt(30))/35;
    pos2 = sqrt(525-70*sqrt(30))/35;
    z = [-pos1; -pos2; pos2; pos1];
    w1 = sqrt(30)*(3*sqrt(30)-5)/180;
    w2 = sqrt(30)*(3*sqrt(30)+5)/180;
    w = [w1   w2   w2   w1];
elseif ngaus == 5
    r70 = sqrt(70);
    pos1 = sqrt(245+14*r70)/21;
    pos2 = sqrt(245-14*r70)/21;
    z = [-pos1; - pos2; 0; pos2; pos1];
    w1 = (7+5*r70)*3*r70/(100*(35+2*r70));
    w2 = -(-7+5*r70)*3*r70/(100*(-35+2*r70));
    w0 = 128/225;
    w = [w1,w2,w0,w2,w1];
elseif ngaus == 6
    z = [0.23861918608319690863050172168066;     0.66120938646626451366139959501991;
         0.93246951420315202781230155449406;    -0.23861918608319690863050172168066;
        -0.66120938646626451366139959501991;    -0.93246951420315202781230155449406];
    w = [0.46791393457269104738987034398801      0.36076157304813860756983351383812, ...
         0.17132449237917034504029614217260      0.46791393457269104738987034398891, ...
         0.36076157304813860756983351383816      0.17132449237917034504029614217271 ];
elseif ngaus == 7
    z = [0.
         0.40584515137739716690660641207692;     0.74153118559939443986386477328078
         0.94910791234275852452618968404784;    -0.40584515137739716690660641207692
        -0.74153118559939443986386477328078;    -0.94910791234275852452618968404784 ];
    w = [0.4179591836734693877551020408166, ...
         0.38183005050511894495036977548841      0.27970539148927666790146777142377, ...
         0.12948496616886969327061143267904      0.38183005050511894495036977548964, ...
         0.27970539148927666790146777142336      0.12948496616886969327061143267912 ];
elseif ngaus == 8
    z = [0.18343464249564980493947614236027;     0.52553240991632898581773904918921
         0.79666647741362673959155393647586;     0.96028985649753623168356086856950
        -0.18343464249564980493947614236027;    -0.52553240991632898581773904918921
        -0.79666647741362673959155393647586;    -0.96028985649753623168356086856950];
    w = [0.3626837833783619829651504492780       0.31370664587788728733796220198797, ...
         0.22238103445337447054435599442573      0.10122853629037625915253135431028, ...
         0.3626837833783619829651504492834       0.31370664587788728733796220198807, ...
         0.22238103445337447054435599442632      0.10122853629037625915253135431015];
elseif ngaus == 9
    z = [0.
         .32425342340380892903853801464336;     .61337143270059039730870203934149
         .83603110732663579429942978806972;     .96816023950762608983557620290365
        -.32425342340380892903853801464336;    -.61337143270059039730870203934149
        -.83603110732663579429942978806972;    -.96816023950762608983557620290365];
    w = [.3302393550012597631645250692903;
         .3123470770400028400686304065887;      .2606106964029354623187428694188
         .18064816069485740405847203124168;    0.81274388361574411971892158110806e-1
        .3123470770400028400686304065836;       .2606106964029354623187428694150
        .18064816069485740405847203124263;     0.81274388361574411971892158110938e-1 ]';
elseif ngaus == 10
    z = [.14887433898163121088482600112972;      .43339539412924719079926594316579
         .67940956829902440623432736511487;      .86506336668898451073209668842349
         .97390652851717172007796401208445;     -.14887433898163121088482600112972
        -.43339539412924719079926594316579;     -.67940956829902440623432736511487
        -.86506336668898451073209668842349;     -.97390652851717172007796401208445];
    w = [.2955242247147528701738929946601;       .26926671930999635509122692156867
         .21908636251598204399553493422796;      .14945134915058059314577633966048
        0.66671344308688137593568809894898e-1;   .2955242247147528701738929946484
        .26926671930999635509122692157323;       .21908636251598204399553493422877
        .14945134915058059314577633965578;      0.66671344308688137593568809893298e-1]';
else
    error('unavailable quadrature')
end
