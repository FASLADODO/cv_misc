function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr, normals, sigmaGM)
%
% function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr,
%                                  normals, sigmaGM)
%
%  minimize sum_i  rho( a'x_i + b + x_i'x_i, sigma)
%  w.r.t. a and b, where a = -2x_0 and b = x_0'x_0 - r^2
% Input:
%  pts: an Nx2 matrix containing 2D x_i, i=1...N
%  initx0: 2-vector initial guess for centre
%  initr: initial guess for radius 
%  normals: an Nx2 matrix of corresponding 2D edge normal directions nrml_i  
% Output
%  x0 - 2x1 fitted circle's center position
%  r  - fitted circle's radius
%  w  - N x 1 robust weights for edgels.
%  maxW  - maximum possible robust weight (i.e. for zero error)
%          This may be used by the calling code to scale the weights 
%          to be between 0 and 1.

% initx0 = circle(1:2);
% initr = circle(3);
% pts = p;

x0 = initx0(:)';  % Make it a row vector.
r = initr;

x = pts(:,1);
y = pts(:,2);

obj_pre = 0;
max_iter = 100;

for i = 1:max_iter
    e = pts-repmat(x0,size(pts,1),1);
    e = e(:,1).^2+e(:,2).^2-r^2;
    obj = sum((e.^2)./(sigmaGM^2+e.^2));
    w = (2*sigmaGM^2)./((sigmaGM^2+e.^2).^2);

    A = [sum(w.*x.^2),sum(w.*x.*y),w'*x;
         sum(w.*x.*y),sum(w.*y.^2),w'*y;
         w'*x,        w'*y,        sum(w)];
    M = [-sum(w.*x.*(x.^2+y.^2));
         -sum(w.*y.*(x.^2+y.^2));
         -sum(w.*(x.^2+y.^2))];
    N = A\M;
    
    x0 = -N(1:2)'/2;
    r = sqrt(x0*x0'-N(3));
    
    if norm(obj-obj_pre) < 1e-6
        break
    else
        obj_pre = obj;
    end    
end

x0 = x0';
maxW = max(w);


% FINISH THIS CODE
