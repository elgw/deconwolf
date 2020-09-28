function getfwhm()

addpath('../')
close all

files(1).name = 'Bead.tif';
files(1).tag = 'Original';
files(end+1).name = 'Final Display of RL.tif';
files(end).tag = 'deLa2';
files(end+1).name = 'Final Display of RL_80it.tif';
files(end).tag = 'deLa2 @80';
files(end+1).name = 'Final Display of RL_160it.tif';
files(end).tag = 'deLa2 @160';
files(end+1).name = 'dw40_Bead.tif';
files(end).tag = 'dw40';
files(end+1).name = 'dw40.bq1_Bead.tif';
files(end).tag = 'dw40 --bq 1';
files(end+1).name = 'dw0_Bead.tif';
files(end).tag = 'dw40 --bq 0';
%files(end+1).name = 'bwX_Bead.tif';
%files(end).tag = 'bwX from newpsf';

% Spatial resolution	deltar: 64.5 nm
% Axial resolution	deltaz: 160 nm
legendstrs = {};
for kk = 1:numel(files)
    file = files(kk).name;
    fprintf('File: %s\n', file);
    I = double(df_readTif(file));
    
%     px = sum(sum(I, 3), 2);
%     py = squeeze(sum(sum(I, 1), 3));
%     pz = squeeze(sum(sum(I, 1), 2));
    %figure, plot(px), hold on, plot(py), plot(pz)
    
    midx = 125.0;
    midy = 128.5;
    midz = 101.75;
    
    % Extract profiles
    d = 45;
    XF = 17; % oversampling
    rx = linspace(midx-d, midx+d, (2*d+1)*XF);
    ry = linspace(midy-d, midy+d, (2*d+1)*XF);
    rz = linspace(midz-d, midz+d, (2*d+1)*XF);
    r1 = ones(size(rx));
    
    px = interpn(I, rx,       midy*r1, midz*r1, 'cubic');
    py = interpn(I, midx*r1,  ry,      midz*r1, 'cubic');
    pz = interpn(I, midx*r1,  midy*r1, rz, 'cubic');
    
    mid = interpn(I, midx, midy, midz)/max(I(:));
    rc = 1-mid;
    files(kk).relContrastPerc = 100*rc;
    %fprintf('RC: %f\n', rc);
    
    levels = [.25, .5, .75];
    imax = max(I(:)); % Use max from whole image
    imax = []; % use max from profile
    for level = levels
        fwhm_xy = (pfwhm(px, level, imax)*64.5/XF + pfwhm(py, level, imax)*64.5/XF)/2;
        fwhm_z = pfwhm(pz, level, imax)*160/XF;
        files(kk).(sprintf('fwhm_xy_%d', 100*level)) = fwhm_xy;
        files(kk).(sprintf('fwhm_z_%d', 100*level)) = fwhm_z;        
        files(kk).mem_Mb = dwGetMem(file)/1000;
        files(kk).time_s = dwGetTime(file);
    end
    
    tab = struct2table(files);
    fprintf('Writing results to tab.csv\n');
    
    df_writeTable(tab, 'tab.csv');
    
    subplot(1,3,1)
    hold on
    plot(rx-midx, px/max(px))
    title('X')
    hold on
    subplot(1,3,2)
    hold on
    plot(ry-midy, py/max(py))
    title('Y')
    subplot(1,3,3)
    hold on
    plot(rz-midz, pz/max(pz))
    title('Z')
    legendstrs{end+1} = files(kk).tag;
    %volumeSlide(I)
    %keyboard
end



subplot(1,3,1)
%legend({'original', 'dwl2', 'deconwolf', 'dw 0', 'dw 1'}, 'location', 'southWest')
legend(legendstrs, 'location', 'southWest')
dprintpdf('radial_profiles', 'publish', 'w', 30, 'h', 6)
end

function fwhm = pfwhm(px, level, imax)
px = px-min(px(:));
if numel(imax > 0)
    px = px/imax;
else
    px = px/max(px);
end

maxx = max(px);
pix = find(px>=level*maxx);
left = pix(1);
right = pix(end);
fwhm = (right-left);
end