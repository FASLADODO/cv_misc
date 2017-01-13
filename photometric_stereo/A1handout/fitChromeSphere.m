function [L] = fitChromeSphere(chromeDir, nDir, chatty)
% [L] = fitChromeSphere(chromeDir, nDir, chatty)
% Input:
%  chromeDir (string) -- directory containing chrome images.
%  nDir -- number of different light source images.
%  chatty -- true to show results images. 
% Return:
%  L is a 3 x nDir image of light source directions.

% Since we are looking down the z-axis, the direction
% of the light source from the surface should have
% a negative z-component, i.e., the light sources
% are behind the camera.
if ~exist('chatty', 'var')
    chatty = false;
end

mask = ppmRead([chromeDir, 'chrome.mask.ppm']);
mask = mask(:,:,1)/255.0;

[chrome_bd_y,chrome_bd_x] = find(mask>0.8);
% chrome center on the image plane (y points down, x points right)
chrome_c = mean([chrome_bd_x,chrome_bd_y]);

chrome_r = chrome_c-[min(chrome_bd_x),min(chrome_bd_y)];
chrome_r = mean(chrome_r);

L = zeros(3,nDir);
reflect_dir = [0,0,-1];

for n = 1:nDir
    fname = [chromeDir,'chrome.',num2str(n-1),'.ppm'];
    im = ppmRead(fname);
    imData(:,:,1) = im(:,:,1)/255.0; %figure(n);imshow(imData)
    [hl_bd_y,hl_bd_x] = find(imData>0.8);

    % find the centroid of the highlighted area
    hl_c = mean([hl_bd_x,hl_bd_y]);
    chrome_surf_n_xy = (hl_c-chrome_c)/chrome_r;
    chrome_surf_norm = [chrome_surf_n_xy,-sqrt(1-sum(chrome_surf_n_xy.^2))];
    % Specular reflectance model: L = 2(N*R)N-R
    L_n = 2*dot(chrome_surf_norm,reflect_dir)*chrome_surf_norm-reflect_dir;
    L(:,n) = L_n';
end

return;

