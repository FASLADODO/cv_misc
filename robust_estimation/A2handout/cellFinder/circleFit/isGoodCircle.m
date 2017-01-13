function [goodCircle] = isGoodCircle(x0, r, w, circleEstimates, nFound)
  % [goodCircle] = isGoodCircle(x0, r, w, normals, ...
  %                                  circleEstimates, nFound)
  % Decide if the circle with parameters x0 and r, with fitted
  % weights w, is to be added to the current set of circles.
  % Input:
  %  x0 2-vector for circle center
  %  r  radius of circle
  %  w  robust weights of edgels supporting this circle,
  %     weights have been scaled to have max possible value 1.
  %  circleEstimates C x 3 array of circle parameters, the first
  %     nFound rows of which specifies the [x0(1), x0(2) r] parameters
  %     for a previously fitted circle for this data set.
  %  nFound the number of previously fitted circles stored in circleEstimates
  % Output:
  %  goodCircle boolean, true if you wish to add this circle to the
  %             list of previously fitted circles.

x0 = x0(:)';  % Decide, row or column
circleEstimates = circleEstimates';

e = 0;
if nFound ~= 0
    e = circleEstimates(1:nFound,1:2)-repmat(x0,nFound,1);
    e = e(:,1).^2+e(:,2).^2-r^2;
    e = sum(e<0);
end

score = sum(floor(w*4)>2)*100/r;
if score > 97 && e <= 0
    goodCircle = 1;
else
    goodCircle = 0;
end

% YOU FINISH THIS