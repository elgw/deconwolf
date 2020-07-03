clear all
close all

s.outFolder = 'set3/';
s.bg = 1000; % Constant added to synthethic image
s.dNoise = 10; % Detector noise power
s.dotSigma = 0.8; % Gaussian size of dots (sigma)
s.nDots = 10000;

mkdir(s.outFolder)

% Target size:
M = 256; N = 256; P = 40;
fprintf('xy-pixels/dot: %.2f\n', M*N/s.nDots);
% PSF size:
pM = 119; pN=119; pP = 2*P-1;
PSF0 = df_readTif('PSF RW.tif');
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

volumeSlide([S(1:M, 1:N, 1:P)/max(S(:)), V(1:M, 1:N, 1:P)/max(V(:))])
