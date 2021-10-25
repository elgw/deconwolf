I = df_readTif('tmr_009.tiff');
outdir = 'crops/';
mkdir(outdir)

cropvalues = [0 200:50:900];

for kk = cropvalues
    J = I(kk+1:size(I)-kk, kk+1:size(I)-kk, :);
    outname = sprintf('%stmr_%03d.tif', outdir, kk);
    if ~isfile(outname)
        df_writeTif(J, outname);
    end
end

mse_dw = zeros(numel(cropvalues), 1);
mse_rl = zeros(numel(cropvalues), 1);
filename0 = sprintf('%sdw_tmr_%03d.tif', outdir, 0);
I0 = df_readTif(filename0);
for ii = 2:numel(cropvalues)
    kk = cropvalues(ii);
    filename = sprintf('%sdw_tmr_%03d.tif', outdir, kk);
    I = df_readTif(filename);
    I0C = I0(kk+1:size(I0)-kk, kk+1:size(I0)-kk, :);
    delta = I - I0C;    
    mse_dw(ii) = mean((delta(:)).^2);
end

T = table(cropvalues(:), mse_dw(:), mse_rl(:) );
T.Properties.VariableNames = {'Cropping', 'MSE_dw', 'MSE_DL2'};

D8 = df_readTif('crops/DL2/RL_800.tif'); D8 = double(D8);
D9 = df_readTif('crops/DL2/RL_900.tif'); D9 = double(D9);
D8c = D8(101:end-100, 101:end-100, :);
deltaD = D8c - D9;


df_writeTif_single(single(abs(deltaD)), 'DL2_900_800.tif');

W8 = df_readTif('crops/dw_tmr_800.tif'); W8 = double(W8);
W9 = df_readTif('crops/dw_tmr_900.tif'); W9 = double(W9);
W8c = W8(101:end-100, 101:end-100, :);
deltaW = W8c - W9;
df_writeTif_single(single(abs(deltaW)), 'DW_900_800.tif');


figure
subplot(2,2,1)
scatter(D8c(:), D9(:), 'k.'), title('Deconvolution Lab2')
corr_dl2 = corr(D8c(:), D9(:));
legend(['\rho = ' sprintf('%.2f', corr_dl2)])
xlabel('Cropped full')
ylabel('crop')
subplot(2,2,2)
scatter(W8c(:), W9(:), 'k.'), title('Deconwolf')
corr_dw = corr(W8c(:), W9(:));
legend(['\rho = ' sprintf('%.2f', corr_dw)])
xlabel('Cropped full')
ylabel('crop')
subplot(2,2,3)
a1 = histogram(deltaD/max(D8c(:)));
set(gca, 'YScale', 'log')
title('N. errors, DL2')
subplot(2,2,4)
a2 = histogram(deltaW/max(W8c(:)));
title('N. errors, DW')
set(gca, 'YScale', 'log')
linkaxes([a1.Parent, a2.Parent])

dprintpdf('subregions_errors', 'driver', {'-dpng', '-dpdf'}, 'w', 20, 'h', 20)

figure
a1 = subplot(2,2,1);
imagesc(D8c(:,:,10))
axis image, colormap gray
title('DL2')
a2 = subplot(2,2,2);
imagesc(W8c(:,:,10))
axis image, colormap gray
title('DW')
a3 = subplot(2,2,3);
imagesc(D9(:,:,10))
axis image, colormap gray
title('DL2 (from crop)')
a4 = subplot(2,2,4);
imagesc(W9(:,:,10))
axis image, colormap gray
title('DW (from crop)')
axes = [a1, a2, a3, a4];
linkaxes(axes)
axis(a1, [0 100, 0, 35])
axis(axes, 'off')

dprintpdf('subregions_zoom', 'driver', {'-dpng', '-dpdf'}, 'w', 25, 'h', 10)

org8 = sprintf('%stmr_%03d.tif', outdir, 800);
org8 = df_readTif(org8);
figure
imagesc(org8(:,:,10)), axis image, colormap gray
hold on
a = 100; b = size(org8,1) - 100;
plot([a, b], [a,a], 'g')
plot([a, b], [b,b], 'g')
plot([a, a], [a,b], 'g')
plot([b, b], [a,b], 'g')
title('RAW')
legend('cropped region')
dprintpdf('subregions_raw',  'driver', {'-dpng', '-dpdf'}, 'w', 15, 'h', 15)