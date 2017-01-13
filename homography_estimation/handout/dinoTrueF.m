%% File: dinoTestF
%% A3 2016 handout code
%% Estimate F matrix from synthetic corresponding points.
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
    cd([userpath_fix,filesep,'csc2503_hw',filesep,'matlabVisTools',filesep]);   %% CHANGE THIS to your startup directory
    run('startup.m');
    cd(dir);
end

reconRoot = [userpath_fix,filesep,'csc2503_hw',filesep,'homography_estimation',filesep,'handout'];  %% CHANGE THIS to the directory you installed A4 handout/
addpath([reconRoot,filesep,'utils']);
cd(reconRoot);

% Random number generator seed:
seed = round(sum(1000*clock));
rand('state',seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1)*1.0e+6);
randn('state',seedn);

nTrial = 10;  %% Number of ransac trials to use


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up cameras
%% The cameras are automatically rotated by projectDino to fixate
%% on the mean of the 3D points.  We do not always want to allow
%% the cameras to move in this fashion (eg. when we change sclZ).
%% So we will compute the rotations for the left and right cameras
%% once and for all, and then use these.
f = 100; % focal length
dLeft = [-50,0,-150]';  % Center of projection for left camera
dRight = [50,0,-150]';  % Center of projection for right camera
%% Compute camera rotations to fixate on Dino's center.
[pLeft,polys,MintLeft,MextLeft] = projectDino(f,dLeft,[],1.0);
Rleft = MextLeft(:,1:3);
[pRight,polys,MintRight,MextRight] = projectDino(f,dRight,[],1.0);
Rright = MextRight(:,1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate data...
sclZ = 1;
%% Dino left image data
[pLeft,polys,MintLeft,MextLeft] = projectDino(f,dLeft,Rleft,sclZ);

%% Dino right image data
[pRight,polys,MintRight,MextRight] = projectDino(f,dRight,Rright,sclZ);

%% Q1
% Left and right cameras projection matrices
P_2 = MintLeft*MextLeft;
P_1 = MintRight*MextRight;

% vector m as the null space of P1
m = null(P_1);
e = P_2*m;
% cross-product matrix
e_x = [    0, -e(3),  e(2);
        e(3),     0, -e(1);
       -e(2),  e(1),     0];
F_0 = e_x*(P_2*pinv(P_1));

%% Show left and right images
hFig = figure(1); clf;
plot(pLeft(1,:), pLeft(2,:), '.b');
axis xy; axis equal;
xlabel('X'); ylabel('Y');
axis([-150,150,-100,100]);
title('Left image of Dino');
pause(0.1);

hFig = figure(2); clf;
plot(pRight(1,:),pRight(2,:),'.b');
axis xy; axis equal;
axis([-150,150,-100,100]);
xlabel('X'); ylabel('Y');
title('Right image of Dino');
pause(0.1);


%% Build correspondence data
clear imPts;
imPts = cat(3,pLeft,pRight);
nPts = size(imPts,2);
if size(imPts,1) == 2
    imPts = cat(1,imPts,ones(1,nPts,2));
end

%% RANSAC for F
seeds = {};
sigma = 2; rho = 2;
for kTrial = 1: nTrial
    %% Test out F matrix on a random sample of 8 points
    idTest = randperm(nPts);
    nTest = min(8,nPts);
    idTest = idTest(1:nTest);
    
    %% Solve for F matrix on the random sample
    [F,Sa,Sf] = linEstF(imPts(:,idTest,1),imPts(:,idTest,2),1);
    
    %% Compute perpendicular error of all points to epipolar lines
    perpErrL = zeros(1,nPts);
    for k = 1:nPts
        lk = imPts(:,k,2)'*F';
        perpErrL(k) = (lk*imPts(:,k,1))/norm(lk(1:2));
    end
    
    %% Detect inliers
    idInlier = abs(perpErrL) < rho*sigma;
    
    %% Count inliers
    nInlier = sum(idInlier);
    if nInlier > 20
        %% Store sets of sampled points with at least 20 inliers
        seed = struct;
        seed.id = idTest;
        seed.idInlier = idInlier;
        seed.nInlier = nInlier;
        seed.F = F;
        
        kSeed = length(seeds)+1
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

%%  Refine estimate of F using all inliers.
F = seeds{ks}.F;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;
sum(idInlier)
%% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1: 10
    %% Fit F using all current inliers
    [F,Sa,Sf] = linEstF(imPts(:,idInlier,1),imPts(:,idInlier,2),1);
    
    %% Compute perpendicular error to epipolar lines
    perpErrL = zeros(1,nPts);
    for k = 1:nPts
        lk = imPts(:,k,2)'*F';
        perpErrL(k) = (lk*imPts(:,k,1))/norm(lk(1:2));
    end
    idInlier = abs(perpErrL) < rho*sigma;
    nInlier = sum(idInlier)
    
    %% If we have the same set of inliers as the previous iteration then stop.
    if all(idInlier == idInlierOld)
        break;
    end
    idInlierOld = idInlier;
end

%% Q2 use pts in the left image
n_col = 16;   % Number of different colours to use.
col = hsv(n_col);  % Colour map.

% cropbox of interest
crop_box_oi = [-50,-50,50,50];

% Create a regularly spaced grid of points
n_space = 10;
grid_lin = linspace(crop_box_oi(1),crop_box_oi(3),n_space);
[grid_x,grid_y] = meshgrid(grid_lin,grid_lin);
% Convert to homogeneous coordinates
grid_pts = [grid_x(:)';
            grid_y(:)';
            ones(size(grid_x(:)'))];

hFig = figure(5);
clf; hold on;
axis([-150,150,-100,100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1),ax(3),ax(2),ax(4)];
title('Epipolar Lines of Grid Points');

d_max = zeros(length(grid_pts),1);

for k = 1:size(grid_pts,2)
    % Plot interest point location corresponding to epipolar line as a "o"
    % in the same colour as the epipolar line.
    plot(grid_pts(1,k),grid_pts(2,k),'o','Color',col(mod(k,n_col)+1,:));
    
    % Plot cropped epipolar line using F0 in the left image using points of
    % interest in the right image.
    lk_0 = grid_pts(:,k)'*F_0';
    epk_0 = cropLineInBox(lk_0(1:2),lk_0(3),crop_box_oi);
    set(line(epk_0(:,1),epk_0(:,2)),'Color',col(mod(k,n_col)+1,:));
    
    % get uncropped epipolar lines using F1
    lk = grid_pts(:,k)'*F';
    epk = cropLineInBox(lk(1:2),lk(3),cropBox);
    set(line(epk(:,1),epk(:,2)),'Color','k','LineStyle',':');
    
    % Convert the endpoints of cropped epipolar line to homogeneous
    % coordinates by adding a column of 1
    epk_h = [epk_0,ones(size(epk_0,1),1)];
    
    % Compute distance from cropped to uncropped epiploar lines
    d = abs(lk*epk_h')./norm(lk(1:2));
    d_max(k) = max(d);
end
max_err = max(d_max)

%% Q3-6 (adjust parameters)
sigma_n = linspace(0,10,31);
err_median = zeros(size(sigma_n));
% Number of samples to get the median of error
n_err = 5;

for n = 1:length(sigma_n)
    % Add Gaussian random noise to the original corresponding points
    rnd_noise = zeros(size(imPts));
    rnd_noise(1:2,:,:) = normrnd(0,sigma_n(n),size(imPts(1:2,:,:)));
    noise_samples = imPts+rnd_noise;

    % Get error samples for each sigma_n
    err_samples = zeros(n_err,1);
    for i = 1:n_err
        % Permutate the total points to get random samples with a size of
        % nPts/err (e.g. 315/5=63), and get these points' IDs
        lines_id = randperm(nPts,8);
        % Get estimate of F from these samples
        [F,~,~] = linEstF(noise_samples(:,lines_id,1),noise_samples(:,lines_id,2),1);
        
        d_max = zeros(size(lines_id));
        for j = 1:size(lines_id,2)
            k = lines_id(j);
            lk = noise_samples(:,k,2)'*F';

            lk_0 = noise_samples(:,k,2)'*F_0';
            epk_0 = cropLineInBox(lk_0(1:2),lk_0(3),crop_box_oi);
            epk_h = [epk_0,ones(size(epk_0,1),1)];

            d = abs(lk*epk_h')./norm(lk(1:2));
            d_max(j) = max(d);
        end
        max_err = max(d_max);
        err_samples(i) = max_err;
    end
    err_median(n) = median(err_samples);
end
figure
plot(sigma_n,err_median);

% rate = diff(err_median)./diff(sigma_n);
% for i = 1:length(sigma_n)-1
%     text(sigma_n(i+1),err_median(i+1),num2str(rate(i)));
% end
xlabel('Sigma_n');
ylabel('Error median');
title('Error median vs. Sigma_n');
%%
