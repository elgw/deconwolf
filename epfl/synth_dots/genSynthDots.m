function genSynthDots(settings)
% Purpose: 
% Generate synthetic image with increasing dot density.
% See how densely they can be packed while still beeing separatable after
% deconvolution.
% 

s.bg = 1000; % Constant added to synthethic image
s.dNoise = 10; % Detector noise power
s.dotSigma = 0.8; % Gaussian size of dots (sigma)
s.nDots = 10000;
s.psfFile = 'PSF RW.tif'; 
s.pixelSizeNM = [130 130 130];
s.imsize = [256 256 40];
s.viewVolume = 0; % View created volume at end

volume = prod(s.imsize.*s.pixelSizeNM/1000) % um3
nucVol = 374; % um3 https://bionumbers.hms.harvard.edu/bionumber.aspx?id=101402&ver=14
nucVol/volume

if exist('settings', 'var')
    s = df_structput(s, settings);
end

s.outFolder = sprintf('sdots_%d/', s.nDots);
s.nDots = round(s.nDots);


mkdir(s.outFolder)
s.logFile = [s.outFolder, 'settings.txt'];
log = fopen(s.logFile, 'a');
df_printstruct(s, log);

% Target size:
M = s.imsize(1); N = s.imsize(2); P = s.imsize(3);
fprintf('xy-pixels/dot: %.2f\n', M*N/s.nDots);
fprintf(log, 'xy-pixels/dot: %.2f\n', M*N/s.nDots);

% PSF size:
pM = 119; pN=119; pP = 2*P-1;
PSF0 = df_readTif(s.psfFile);
while size(PSF0,3) > pP
    PSF0 = PSF0(2:end-1, 2:end-1, 2:end-1);
end
PSF0 = PSF0/sum(PSF0(:));

df_writeTif_single(PSF0, [s.outFolder 'psf.tif']);

T = zeros(M+(pM+1)/2, N+(pN+1)/2, P+(pP+1)/2); 
PSF = zeros(size(T));
PSF(1:size(PSF0,1), 1:size(PSF0,2), 1:size(PSF0,3)) = PSF0;
for kk = 1:3
    PSF = circshift(PSF, -(size(PSF0, kk)-1)/2, kk);
end


%% Positions
X = [ 1 + rand(s.nDots, 1)*(M-1), ...
      1 + rand(s.nDots, 1)*(N-1), ...
      1 + rand(s.nDots, 1)*(P-1)];

out = fopen([s.outFolder 'coords.csv'], 'w');
fprintf(out, 'X, Y, Z\n');
for kk = 1:size(X,1)
    fprintf(out, '%f, %f, %f\n', X(kk,1), X(kk,2), X(kk,3));
end
fclose(out);

save([s.outFolder 'X.mat'], 'X', 's');
copyfile('genSynthDots.m', s.outFolder);
copyfile('runDw', s.outFolder);

%% Add intensity and sigma/size
% X: x, y, z, intensity, sigma_x, sigma_y, sigma_z
X = [X, 10000+rand(s.nDots,1)*10000, s.dotSigma*ones(s.nDots,3)];
V = df_blit3(T, [], X', 2);
df_writeTif_single(single(V(1:M, 1:N, 1:P)), [s.outFolder 'Ground Truth.tif']);
PSF = PSF/max(PSF(:)); % To not change total number of photons
S = ifftn(fftn(V).*fftn(PSF)); % Convolve

S = S + s.bg; % Background
S = S + poissrnd(S); % Poisson noise
S = S + s.dNoise*randn(size(S)); % Detector noise
% Write as 32-bit float
df_writeTif_single(single(S(1:M, 1:N, 1:P)), [s.outFolder 'sdots.tif']);

if s.viewVolume
    volumeSlide([S(1:M, 1:N, 1:P)/max(S(:)), V(1:M, 1:N, 1:P)/max(V(:))])
end
fclose(log);
