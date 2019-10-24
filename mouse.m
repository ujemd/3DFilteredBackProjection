clear all
close all
clc
%% read data
fname = 'mouse.tif';
info = imfinfo(fname);
num_images = numel(info);
Nphi = num_images;

%% read projections
projs = zeros(info(1).Height, info(1).Width,num_images);
for idx = 1:num_images
    projs(:,:,idx) = im2double(imread(fname,idx));
end

%% visualize projections
for idx = 1:num_images
    imshow(projs(:,:,idx),[])
    pause(0.001)
end

%% get sinograms
Nprojs = size(projs,1);
sinograms = zeros(size(projs,2),Nphi,Nprojs);
for idx = 1:Nprojs
    sinograms(:,:,idx) = projs(idx,:,:);
end

%% visualize sinograms
for idx = 1:size(projs,1)
    imagesc(sinograms(:,:,idx))
    axis xy; axis image; colorbar; colormap gray;
    pause(0.001)
end

%% Geometry setup
Nx = 270;
Ny = 104;
Nr = 290;
Nphi = 180;

% Projection dimensions
xLen = 54010; %microns
yLen = 20800; %microns
phiLen   = pi;
rLen =   58011; % length of detector

phiVec   = (0	:   Nphi   - 1)   *    phiLen   / Nphi;
rVec =   (0.5	:   Nr -   0.5)   *    rLen /   Nr - rLen / 2;
xVec =   (0.5	:   Nx -   0.5)   *    xLen /   Nx - xLen / 2;
yVec =   (0.5	:   Ny -   0.5)   *    yLen /   Ny - yLen / 2;

%% Filter Sinograms
q = zeros(size(sinograms));
for idx = 1:Nprojs
    q(:,:,idx) = rampfilter(sinograms(:,:,idx), rVec ,'fourier1');
end

%% visualize sinograms
for idx = 1:size(projs,1)
    subplot(1,2,1); imagesc(sinograms(:,:,idx)); title('Sinogram');
    axis xy; axis image; colorbar; colormap gray;
    subplot(1,2,2); imagesc(q(:,:,idx)); title('Filtered');
    axis xy; axis image; colorbar; colormap gray;
    pause(0.001)
end

%% regular backprojection

bp = zeros(Ny,Nx,Nprojs);
for idx = 1:Nprojs
    bp(:,:,idx) = backproject(sinograms(:,:,idx), rVec , phiVec , xVec , yVec ,'linear');
end

%% filtered backprojection

fbp = zeros(Ny,Nx,Nprojs);
for idx = 1:Nprojs
    fbp(:,:,idx) = backproject(q(:,:,idx), rVec , phiVec , xVec , yVec ,'linear');
end

%% iradon

%base = zeros(Ny,Nx,Nprojs);
for idx = 1:Nprojs
    base(:,:,idx) = iradon(sinograms(:,:,idx), 0:Nphi-1,'linear','Ram-Lak',1,Nx);
end

%% visualize

%% transverse plane regular backprojection
for idx = 1:Nprojs
    imagesc(bp(:,:,idx)); axis xy; axis image; colorbar; colormap gray;
    pause(0.001);
end

%% coronal plane regular bakcprojection
for idx = 1:size(bp,1)
    imagesc(reshape(bp(idx,:,:),[Nx Nprojs])); axis xy; axis image; colorbar; colormap gray;
    pause(0.001);
end

%% transverse plane filtered backprojection
for idx = 1:Nprojs
    imagesc(abs(fbp(:,:,idx))); axis xy; axis image; colorbar; colormap gray;
    pause(0.001);
end

%% coronal plane filtered backprojection
for idx = 1:size(fbp,1)
    imagesc(abs(reshape(fbp(idx,:,:),[Nx Nprojs]))); axis xy; axis image; colorbar; colormap gray;
    M(idx) = getframe(gcf);
    pause(0.01);
end

%% transverse iradon
for idx = 1:size(base,1)
    imagesc(abs(squeeze(base(idx,:,:)))); axis xy; axis image; colorbar; colormap gray;
    pause(0.001);
end

%% create movie
movie(figure, M)
%% save movie
v = VideoWriter('mouse.avi');
open(v);
writeVideo(v,M);
close(v);

%% save images
fbp2 = abs(fbp);
fbp2 = fbp2/max(fbp2(:));
outputFileName = 'result.tif';
for idx=1:size(fbp2,1)
   imwrite(reshape(fbp2(idx,:,:),[Nx Nprojs]), outputFileName, 'WriteMode', 'append', 'Compression','none');
end