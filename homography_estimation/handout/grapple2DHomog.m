%% File: grappleFmatrix
%% A3 2016 handout code
%% Uses RANSAC to estimate H matrix from corresponding points.
%%
%% ADJ

clear
close all
FALSE = 1 == 0;
TRUE = ~FALSE;
userpath_fix = userpath;
userpath_fix = userpath_fix(1:end-1);

global matlabVisRoot

%% We need to ensure the path is set for the iseToolbox.
if ~isempty(length(matlabVisRoot))
    dir = pwd;
    cd([userpath_fix,filesep,'csc2503_hw',filesep,'matlabVisTools',filesep]);
    run('startup.m');
    cd(dir);
end

reconRoot = [userpath_fix,filesep,'csc2503_hw',filesep,'homography_estimation',filesep,'handout'];
addpath([reconRoot,filesep,'data',filesep,'wadham']);
addpath([reconRoot,filesep,'utils']);
cd(reconRoot);

% Random number generator seed:
seed = round(sum(1000*clock));
rand('state',seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1)*1.0e+6);
randn('state',seedn);

nTrial = 10;  % Number of ransac trials to try.

%% Wadham left image: use  wadham/001.jpg
imPath = ['data',filesep,'wadham',filesep]; fnameLeft = '001';
im = imread([imPath,fnameLeft],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imLeft = imDwn;

%% Read correspondence data
load(['data',filesep,'wadham',filesep,'corrPnts5'])
%% Wadham right image data/wadham002-5.jpg use for corrPnts2-5 respectively
fnameRight = '005';
im = imread([imPath,fnameRight],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imRight = imDwn;

clear imPts;
imPts = cat(3,im_pos1',im_pos2');
nPts = size(imPts,2);
if size(imPts,1)==2
    imPts = cat(1,imPts,ones(1,nPts,2));
end

%% RANSAC for H
seeds = {};
sigma = 2.0; rho = 2;
for kTrial = 1:nTrial
    %% Test out H matrix on a random sample of 4 points
    idTest = randperm(nPts);
    nTest = min(4,nPts);
    idTest = idTest(1:nTest);
    
    %% Solve for H matrix on the random sample
    [H,Sa] = linEstH(imPts(:,idTest,1),imPts(:,idTest,2),1);

    %% Compute error between 2 points
    dist_err = zeros(1,nPts);
    for k = 1:nPts
        p = imPts(:,k,1);
        q = imPts(:,k,2);
        qhat = H*p;
        phat = pinv(H)*q;
        qhat = qhat/qhat(end);
        phat = phat/phat(end);
        %cross(imPts(:,k,2),q_k)
        dist_err(k) = norm(q-qhat)+norm(p-phat);
    end
    
    %% Detect inliers
    idInlier = abs(dist_err)<rho*sigma;
    
    %% Count inliers
    nInlier = sum(idInlier);
    if nInlier > 20
        %% Store sets of sampled points with at least 20 inliers
        seed = struct;
        seed.id = idTest;
        seed.idInlier = idInlier;
        seed.nInlier = nInlier;
        seed.H = H;
        
        kSeed = length(seeds)+1;
        seeds{kSeed} = seed;
    end
end
%% Done RANSAC trials

%% Extract best solution
nInliers = zeros(1,length(seeds));
for ks = 1:length(seeds)
    nInliers(ks) = seeds{ks}.nInlier;
end
[nM,ks] = max(nInliers);
nInliers(ks)

%%  Refine estimate of H using all inliers.
H = seeds{ks}.H;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;
sum(idInlier)
%% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1:10
    %% Fit F using all current inliers
    [H,Sa] = linEstH(imPts(:,idInlier,1),imPts(:,idInlier,2),1);
    
   %% Compute error
    dist_err = zeros(1,nPts);
    for k = 1:nPts
        p_k = imPts(:,k,1);
        q_k = H*p_k;
        q_k = q_k/q_k(end);
        dist_err(k) = norm(imPts(:,k,2)-q_k);
        %cross(imPts(:,k,2),q_k)
    end
    idInlier = abs(dist_err)<rho*sigma;
    nInlier = sum(idInlier)
    
    %% If we have the same set of inliers as the previous iteration then stop.
    if all(idInlier==idInlierOld)
        break;
    end
    idInlierOld = idInlier;
end

%% Show left image
SUPERIMPOSE = TRUE;
hFig = figure(1);
clf;
if SUPERIMPOSE
    image(imLeft);
    colormap(gray(256));
end
%resizeImageFig(hFig,size(imDwn),1);
title(['Left Image',num2str(fnameLeft)])

%% Show warped image
hFig = figure(2);
clf;
if SUPERIMPOSE
    image(homogWarp(imLeft,inv(H),size(imDwn)));
    colormap(gray(256));
end
%resizeImageFig(hFig,size(imDwn),1);
title('Warped Left Image')

%% Show right image
SUPERIMPOSE = TRUE;
hFig = figure(3);
clf;
if SUPERIMPOSE
    image(imRight);
    colormap(gray(256));
end
%resizeImageFig(hFig,size(imDwn),1);
title(['Right Image',num2str(fnameRight)])

