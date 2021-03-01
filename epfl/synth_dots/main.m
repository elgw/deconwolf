function main()
clear all
close all


%% Settings
N = [1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
N = [11000 12000 13000 14000 15000];
N = [16000 17000 18000 19000 20000];
N = (1:25)*1000;
N = round(N);

volume_um3 = prod(1/1000*[130 130 130].*[256 256 40]);
cell_volume_um3 = 3000;
fprintf('Image volume: %f um3, Cell volume: %f um3\n', volume_um3, cell_volume_um3);
for kk = 1:numel(N)
    fprintf('N=%d, dots/cell=%.1f\n', N(kk), N(kk)/volume_um3*cell_volume_um3);
end

s.gen = 0; % Generate new images?
s.genSynth = 0;
s.nIter = 150; % Number of iterations of deconwolf

%% Run
if s.gen == 1
    if s.genSynth == 1
        warning('This will remove the originals!');
        keyboard
        % Syntetic image
        for kk = 1:numel(N)
            s1.nDots = N(kk);
            genSynthDots(s1);
        end
    end
    
    % Deconvolve
    for kk = 1:numel(N)
        folder = sprintf('sdots_%d/', N(kk));
        image = sprintf('%ssdots.tif', folder);
        psf = sprintf('%spsf.tif', folder);
        command = sprintf('LD_LIBRARY_PATH=''''; dw --overwrite --iter %d --prefix %d %s %s', s.nIter, s.nIter, image, psf);
        [status, message] = system(command);
        if status ~= 0
            disp(message)
        end
    end
end

% 20201020 - switched from analyse to analyse2
outfolder =  'csv20200120/';
mkdir(outfolder);
method = @analyse2;
files = {'Ground Truth.tif', 'sdots.tif'}
for niter = [25, 50, 100, 125, 150]
    files{end+1} = sprintf('%d_sdots.tif', niter);
end
for niter = [50, 100, 200]
    files{end+1} = sprintf('dl2_%d_sdots.tif', niter);
end

fprintf('Files:\n');
for kk = 1:numel(files)       
    fprintf('newtab%d : %s\n', kk, files{kk});
end
pause
for kk = 1:numel(files)       
    fprintf('Pattern: %s\n', files{kk});
    tab = method(s, N, files{kk});    
    df_writeTable(tab, sprintf('%s/newtab%d.csv', outfolder, kk));
end
disp('Done')



end

function tab = analyse2(s, N, filename)
%% Analyze the images, for all N
for kk = 1:numel(N)
    folder = sprintf('sdots_%d/', N(kk));
    settings = [];
    settings.folder = folder;
    settings.file = filename; 
    settings.inputfile = sprintf('sdots.tif');
    status = getmse2(settings);
    % Add number of dots
    a = strsplit(status.folder, '_');
    status.N = str2num(a{2}(1:end-1));
    % Fix folder name (to constant width)
    fldr = '             ';
    fldr(end-numel(status.folder)+1:end) = status.folder;
    status.folder = fldr;
    
    if exist('tab', 'var')
        tab = [tab; struct2table(status)];
    else
        tab = struct2table(status);
    end
end
end

function tab = analyse(s, N, filename)
%% Analyze the images
for kk = 1:numel(N)
    folder = sprintf('sdots_%d/', N(kk));
    settings = [];
    settings.folder = folder;
    settings.file = filename;
    settings.inputfile = sprintf('sdots.tif');
    status = getmse(settings);
    % Add number of dots
    a = strsplit(status.folder, '_');
    status.N = str2num(a{2}(1:end-1));
    % Fix folder name (to constant width)
    fldr = '             ';
    fldr(end-numel(status.folder)+1:end) = status.folder;
    status.folder = fldr;
    
    if exist('tab', 'var')
        tab = [tab; struct2table(status)];
    else
        tab = struct2table(status);
    end
end
end