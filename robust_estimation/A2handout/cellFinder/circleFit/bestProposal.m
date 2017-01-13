function [circle] = bestProposal(circles, sigmaGM, normals, p)
% [circle] = bestProposal(circles, sigmaGM, normals, p)
% Chose the best circle from the proposed circles.
% Input
%  circles K x 3 matrix each row containing [x0(1), x0(2), r] for 
%          the center and radius of a proposed circle.
%  sigmaGM - the robust scale parameter 
%  normals - P x 2 edgel normal data
%  p       - P x 2 edgel position data.
% Output
%  circle  1 x 3  best circle params [x0(1), x0(2), r] from the proposals

rho_total = zeros(size(circles,1),1);
for i = 1:size(circles,1)
    e = p-repmat(circles(i,1:2),size(p,1),1);
    e = e(:,1).^2+e(:,2).^2-circles(i,3)^2;
    rho_total(i,:) = sum((e.^2)./(sigmaGM^2+e.^2));
end
[~,id] = min(rho_total);
circle = circles(id,:);

  % YOU COMPLETE THIS
