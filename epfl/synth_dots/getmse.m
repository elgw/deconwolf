function status = getmse(settings)
% function status = getmse(struct settings)
% Required fileds of settings: 'folder', 'file'
% Optional: 'Tag'
s.Tag = 'No Tag';

s = df_structput(s, settings);
status = s;

addpath('../')

% Note about SSIM:
% in sage1701 they probably used https://www.imatest.com/docs/ssim/

% We don't scale the reference image
Iref = df_readTif([s.folder 'Ground Truth.tif']);
Iref = double(Iref);

dwfile = [s.folder s.file];
rfile = [s.folder s.inputfile];
fprintf('file: %s\n', dwfile);

I0 = double(df_readTif(rfile)); % Non-deconvolved file

if isfile(dwfile)
    I = df_readTif(dwfile);
else
    warning('Can''t open %s\n', dwfile);
    I = -1*ones(size(I0));
end

if contains(dwfile, 'dl2')
    warning('Shifting dots based on the image file name');
    I = circshift(I, -1, 1);
    I = circshift(I, -1, 2);
end

I = double(I);
I = I*mean(Iref(:))/mean(I(:)); % For easier opimization


[scaling, mse] = fminunc(@(x) mseval(Iref, I*x), 0.1);
%fprintf('mse(%s,ref) = %f\n', s.Tag, mse);
status.mse = mse;

[scaling, mse] = fminunc(@(x) mseval(Iref, I0*x), 0.1);
%fprintf('mse(%s,ref) = %f\n', s.Tag, mse);
status.mse_raw = mse;

[a, b] = getDots(s, I, Iref);
status.DWnDotsIsLMAX = a*100;
status.DWnDotsIsLMAX2 = b*100;
%keyboard
[a, b] = getDots(s, I0, Iref);

status.nDotsIsLMAX = a*100;
status.nDotsIsLMAX2 = b*100;


status.time_s = dwGetTime(dwfile);
status.mem_kb = dwGetMem(dwfile);
%status.ssim = ssim(I, Iref);
%status.psnr = psnr(I, Iref);

end


function v = mseval(A, B)
v = sqrt(mean((A(:)-B(:)).^2));
end


function [N, N2]  = getDots(s, I, Iref)
% Returns how many of the true dots that are local maximas
%
% TODO: over all thresholds
% TODO: fix potential shifts with Deconvolution Lab2 ... why?

METHOD_NDOTS = 2;

method = METHOD_NDOTS;

D = load([s.folder 'X.mat']);
D = D.X; % Coordinates of the true dots
[x, y,z, p, Map] = getLocalMaximas(I);
[xr, yr,zr, pr, Mapr] = getLocalMaximas(Iref);

if method == METHOD_NDOTS
N = numel(x)/numel(xr);
N2 = numel(x)/numel(xr);
return
end


if 0
    figure,
    plot(x, y, 'ro')
    hold on
    plot(xr+1, yr+1, 'gx')
    figure,
    plot(x, z, 'ro')
    hold on
    plot(xr+1, zr, 'gx')
end



N = sum(Map(:) == 1 & Mapr(:) == 1)/sum(Mapr(:));

Map2 = max(Map, [], 3); Mapr2 = max(Mapr, [], 3);
N2 = sum(Map2(:) == 1 & Mapr2(:) == 1)/sum(Mapr2(:));


return
keyboard

Di = round(D(:,1:3));
trueMap = 0*Map;
trueMap(sub2ind(size(I), Di(:,1), Di(:,2), Di(:,3))) = 1;
trueMap = clearBoarders(trueMap, 4, 0); % Same as in getLocalMaximas

Di = Di(min(Di, [], 2) > 0, :);
whos D
whos Di
idx = sub2ind(size(I), Di(:,1), Di(:,2), Di(:,3));
N = sum(Map(:) == 1 & trueMap(:) == 1);
N = N/sum(trueMap(:));

if 0
    figure
    imagesc(sum(I,3)), title('Image')
    figure, imagesc(sum(Map,3)), title('Maximas')
    Map2 = 0*Map; Map2(idx) = 1;
    figure, imagesc(sum(Map2, 3)); title('True locations')
end

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
Map = (Image > D) & (Image > (mean(Image(:))+std(Image(:))));
Pos = find(Map == 1);
[PX, PY, PZ]=ind2sub(size(Image), Pos);

end