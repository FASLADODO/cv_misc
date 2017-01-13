function [depth] = getDepthFromNormals(n,mask)
% [depth] = getDepthFromNormals(n, mask)
%
% Input:
%    n is an [N, M, 3] matrix of surface normals (or zeros
%      for no available normal).
%    mask logical [N,M] matrix which is true for pixels
%      at which the object is present.
% Output
%    depth an [N,M] matrix providing depths which are
%          orthogonal to the normals n (in the least
%          squares sense).
%
imsize = size(mask);
[N,M] = size(mask);

nx = n(:,:,1);
ny = n(:,:,2);
nz = n(:,:,3);
nx = nx(:);
ny = ny(:);
nz = nz(:);

A = sparse(2*N*M,N*M);
A = spdiags(nz(:),0,A);
A = spdiags([zeros(N,1);-nz(1:end-N)],N,A);
A = spdiags(nz(:),-N*M,A);
A = spdiags([zeros(1,1);-nz(1:end-1)],-(N*M-1),A);

v = sparse([nx(:);ny(:)]);

z = A\v;
ground_idx = find(z~=0,1);
z = z-z(ground_idx);
depth = full(reshape(z,imsize));
depth = depth.*mask;

  % YOU NEED TO COMPLETE THIS.
