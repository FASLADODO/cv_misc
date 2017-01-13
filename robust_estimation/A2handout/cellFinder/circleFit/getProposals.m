function [circles] = getProposals(normals, p, numGuesses)
% [circles] = getProposals(normals, p, numGuesses)
% Attempt to produce up to numGuesses circle proposals from
% the edgel data p and normals.  For typical data sets
% we will be able to produce numGuesses proposals.  However,
% on some datasets, say with only a few edgels, we may not be
% able to generate any proposals.  In this case size(circles,1)
% can be zero.
% Input:
%  normals - N x 2 edgel normals
%  p         N x 2 edgel positions
%  numGuesses - attempt to propose this number of circles.
% Return:
%   circles a P x 3 array, each row contains [x0(1) x0(2) r]
%           with 0 <= P <= numGuesses.

n_r = 10;
r = linspace(6,20,n_r)';
circ_prop = zeros(n_r,2,size(p,1));

for i = 1:size(p,1)
    circ_prop(:,:,i) = repmat(p(i,:),n_r,1)+r*normals(i,:);
end
circ_prop = ipermute(circ_prop,[1,3,2]);
circ_prop = reshape(circ_prop,n_r*size(p,1),2);
circ_prop = [circ_prop,repmat(r,size(p,1),1)];

x_bounds = [min(circ_prop(:,1)),max(circ_prop(:,1))];
y_bounds = [min(circ_prop(:,2)),max(circ_prop(:,2))];
x_step = x_bounds(1):5:x_bounds(2);
y_step = y_bounds(1):5:y_bounds(2);

peaks = zeros((size(y_step,2)-1)*(size(x_step,2)-1),3);
n = 0;
for i = 1:size(x_step,2)-1
    for j = 1:size(y_step,2)-1
        pts = sum(circ_prop(:,1)>x_step(i) & circ_prop(:,1)<x_step(i+1) &...
                  circ_prop(:,2)>y_step(j) & circ_prop(:,2)<y_step(j+1));
        n = n+1;
        peaks(n,:) = [pts,i,j];
    end
end
peaks = sortrows(peaks,-1);
num_prop = min(sum(peaks(:,1)>1),numGuesses);
peaks = peaks(1:num_prop,:);

circles = zeros(num_prop,3);
for i = 1:num_prop
    x_idx = peaks(i,2);
    y_idx = peaks(i,3);
    idx = circ_prop(:,1)>x_step(x_idx) & circ_prop(:,1)<x_step(x_idx+1) & ...
          circ_prop(:,2)>y_step(y_idx) & circ_prop(:,2)<y_step(y_idx+1);
    circles_guesses_sub = circ_prop(idx,:);
    circles(i,:) = mean(circles_guesses_sub);
end


  % YOU NEED TO FILL IN CODE HERE.
  
