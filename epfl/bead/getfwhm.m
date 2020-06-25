
function getfwhm()
close all

files = {'dw_Bead.tif', 'Final Display of RL.tif', 'Bead.tif'}

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
    d = 30;
    px = interpn(I, midx-d:midx+d,        midy*ones(1, 2*d+1), midz*ones(1, 2*d+1));
    py = interpn(I, midx*ones(1, 2*d+1),  midy-d:midy+d,       midz*ones(1, 2*d+1));
    pz = interpn(I, midx*ones(1, 2*d+1),  midy*ones(1, 2*d+1), midz-d:midz+d      );
    
    mid = interpn(I, midx, midy, midz)/max(I(:));
    rc = 1-mid;
    
    fprintf('RC: %f\n', rc);
    
    fwhm_x = pfwhm(px)*64.5;
    fwhm_z = pfwhm(pz)*160;            
    
    fprintf('fwhm radial: %f, axial: %f\n', fwhm_x, fwhm_z);        
    
    if 0
    figure
    imagesc(I(:,:,round(midz))), axis image, colormap gray
    hold on, plot(midy, midx, 'rx')
    keyboard
    end
    
    figure
    plot(-d:d, px)
    hold on     
    plot(-d:d, py)
    plot(-d:d, pz)
    title(file, 'interpreter', 'none')
    legend({'x', 'y', 'z'})
       %volumeSlide(I)
       %keyboard
end

end

function fwhm = pfwhm(px)
    maxx = max(px);
    pix = find(px>=0.5*maxx);
    left = pix(1);
    right = pix(end);
    fwhm = (right-left);
    end