function getmse()
addpath('../')
% Note: in sage1701 they probably used https://www.imatest.com/docs/ssim/
%

%psf = df_readTif('psf.tif');

files(1).Name = 'sdots.tif';
files(end).Tag = 'input: w. noise';
files(end+1).Name = 'dw50_sdots.tif';
files(end).Tag = 'dw50';
files(end+1).Name = 'dw100_sdots.tif';
files(end).Tag = 'dw100';
%files(end+1).Name = 'dw150_sdots.tif';
%files(end).Tag = 'dw150';
files(end+1).Name = 'dw200_sdots.tif';
files(end).Tag = 'dw200';
files(end+1).Name = 'dw250_sdots.tif';
files(end).Tag = 'dw250';
files(end+1).Name = 'Final Display of RL 250 iter.tif';
files(end).Tag = 'DeLa2 250';
files(end+1).Name = 'Ground Truth.tif';
files(end).Tag = 'Ground Truth';

% We don't scale the reference image
Iref = df_readTif('Ground Truth.tif');
Iref = double(Iref); 

for ff = 1:numel(files)    
    file = files(ff).Name;
    fprintf('file: %s\n', file);
    tag = files(ff).Tag;
    I = df_readTif(file);
    I = double(I);    
    I = I*mean(Iref(:))/mean(I(:)); % For easier opimization
    [scaling, mse] = fminunc(@(x) mseval(Iref, I*x), 0.1);
    fprintf('mse(%s,ref) = %f\n', tag, mse);
    files(ff).nDotsIsLMAX = getDots(I)*100;
    %files(ff).mse = mse;
    %files(ff).time_s = dwGetTime(file);
    %files(ff).mem_kb = dwGetMem(file);
    %files(ff).ssim = ssim(I, Iref);
    %files(ff).psnr = psnr(I, Iref);
end

tab = struct2table(files);
% writetable(tab, 'tab.csv')
% save as table circumvent whimsical rounding
df_writeTable(tab, 'tab.csv');
end


function v = mseval(A, B)
 v = sqrt(mean((A(:)-B(:)).^2));
end


function N  = getDots(I)

D = load('X.mat');
D = D.X; % Coordinates of the true dots
[x, y,z, p, Map] = getLocalMaximas(I);

Di = round(D(:,1:3));
Di = Di(min(Di, [], 2) > 0, :);
whos D
whos Di
idx = sub2ind(size(I), Di(:,1), Di(:,2), Di(:,3));
N = sum(Map(idx))/size(Di,1);
 
end


function [PX, PY, PZ,  Pos, Map] = getLocalMaximas(Image)
conn26 = 0;
if conn26
sel = ones(3,3,3);
else
sel = zeros(3,3, 3);
sel(:, 2, 2) = 1;
sel(2, :, 2) = 1;
sel(2, 2, :) = 1;
end
sel(2,2, 2) = 0;

D = imdilate(Image, strel('arbitrary', sel));
D = clearBoarders(D, 4, Inf);
Pos = find(Image>D);
[PX, PY, PZ]=ind2sub(size(Image), Pos);
Map = Image>D;
end