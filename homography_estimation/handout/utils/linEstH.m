function [H,Sa] = linEstH(left, right, NUM_RESCALE)
% [H, Sa] = linEstH(left, right, NUM_RESCALE)
% Estimate the homography matrix H from two 3 x n matrices of
% corresponding points left and right.
% Here left(:,k) x (H * right(:,k)) apprx= 0 where x denotes
% the 3D cross product.
% NUM_RESCALE (default TRUE) uses Hartley's rescaling. Always use
% rescaling, unless you wish to show how badly the un-normalized
% algorithm works.
% Returns H along with the singular values Sa of the 2nPts x 9 homogeneous
% linear system for H.

%%
if nargin < 3
    NUM_RESCALE = 1;
end

nPts = size(left,2);
if nPts < 4 || nPts ~= size(right,2)
    fprintf(2,'lineEstH: Innappropriate number of left and right points.');
    H = [];
    return;
end

if size(left,1) == 2
    left = [left;
            ones(1,nPts)];
else % Normalize to pixel coords
    left = left./repmat(left(3,:),3,1);
end
if size(right,1) == 2
    right = [right;
             ones(1,nPts)];
else % Normalize to pixel coords
    right = right./repmat(right(3,:),3,1);
end

imPts = cat(3,left,right);

%% Rescale image data for numerical stability.
if NUM_RESCALE
    Knum = repmat(eye(3),[1,1,2]);
    %%% Rescale for numerical stability
    mn = sum(imPts(1:2,:,:),2)/nPts;
    mns = reshape(mn,[2,1,2]);
    var = sum(sum((imPts(1:2,:,:)-repmat(mns,[1,nPts,1])).^2,2)/nPts,1);
    %% Scale image points so that sum of variances of x and y = 2.
    scl = sqrt(2./var(:));
    %% Sanity: varScl =  var .* reshape(scl.^2, [1 1 2]) % Should be 2
    
    %% Scale so x and y variance is roughly 1, translate so image mean (x,y) is zero.
    Knum(1:2,3,:) = -mn;
    Knum(1:2,:,:) = Knum(1:2,:,:).*repmat(reshape(scl,[1,1,2]),[2,3,1]);
    for kIm = 1:2
        imPts(:,:,kIm) = reshape(Knum(:,:,kIm),3,3)*imPts(:,:,kIm);
    end
    %% Sanity check
    % sum(imPts(1:2,:,:),2)/nPts  % Should be [0 0]'
    % sum(sum(imPts(1:2,:,:).^2,2)/nPts,1) % Should be 2.
end

%% Make constraint matrix A.
%% You fill in below here.
left = reshape(imPts(:,:,1),[3,nPts]);
right = reshape(imPts(:,:,2),[3,nPts]);

A = -left';
[n,m] = size(A); % e.g. 8x3

% Left part of matrix A
B = permute(reshape(A',m,1,[]),[2,1,3]);
C = [B;zeros(1,m,size(B,3))];
C = [C,flipud(C)];
A_1 = reshape(ipermute(C,[1,3,2]),n*2,m*2);

% Right part of matrix A
D = right(1:2,:);
D = repmat(D(:),1,m);
E = reshape(repmat(left(:)',2,1),2,m,[]);
E = reshape(ipermute(E,[1,3,2]),[],m);
A_2 = D.*E;

% Solve linear least square using svd
A = [A_1,A_2];
[~,Sa,Va] = svd(A);

H = reshape(Va(:,end),[3,3])';
Sa = diag(Sa);

%% Sanity
%cross(right(:,1),H*left(:,1))

%% Undo the renormalization
if NUM_RESCALE
    H = reshape(Knum(:,:,2),3,3)\H*reshape(Knum(:,:,1),3,3);
end

H = H./norm(H);
%sumsqr(H)
