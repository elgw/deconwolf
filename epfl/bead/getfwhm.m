
function getfwhm()
close all

files = {'Bead.tif', 'Final Display of RL.tif', 'dw_Bead.tif', 'dw0_Bead.tif', 'dw1_Bead.tif'}

% Spatial resolution	deltar: 64.5 nm
% Axial resolution	deltaz: 160 nm

for kk = 1:numel(files)    
    file = files{kk};
    fprintf('File: %s\n', file);
    I = double(df_readTif(file));
    
    px = sum(sum(I, 3), 2);
    py = squeeze(sum(sum(I, 1), 3));
    pz = squeeze(sum(sum(I, 1), 2));
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
    
    fprintf('RC: %f\n', rc);
    
    imax = max(I(:)); % Use max from whole image
    imax = []; % use max from profile
    fwhm_xy = (pfwhm(px, imax)*64.5/XF + pfwhm(py, imax)*64.5/XF)/2;
    fwhm_z = pfwhm(pz, imax)*160/XF;            
    
    fprintf('fwhm radial: %f, axial: %f\n', fwhm_xy, fwhm_z);        
    
    if 0
    figure
    imagesc(I(:,:,round(midz))), axis image, colormap gray
    hold on, plot(midy, midx, 'rx')
    keyboard
    end
    
    
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
    
       %volumeSlide(I)
       %keyboard
end

subplot(1,3,1)
legend({'original', 'dwl2', 'deconwolf', 'dw 0', 'dw 1'}, 'location', 'southWest')
dprintpdf('radial_profiles', 'publish', 'w', 30, 'h', 6)
end

function fwhm = pfwhm(px, imax)
px = px-min(px(:));
if numel(imax > 0)
    px = px/imax;
else
    px = px/max(px);
end

    maxx = max(px);
    pix = find(px>=0.5*maxx);
    left = pix(1);
    right = pix(end);
    fwhm = (right-left);
    end