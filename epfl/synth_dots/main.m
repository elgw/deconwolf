function main()
clear all
close all


%% Settings
N = [1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
N = [11000 12000 13000 14000 15000];
N = [16000 17000 18000 19000 20000];
N = (1:25)*1000;
N = round(N);
 
s.gen = 0; % Generate new images?
s.nIter = 50; % Number of iterations of deconwolf

%% Run
if s.gen == 1
    % Syntetic image
    for kk = 1:numel(N)
        s1.nDots = N(kk);
        genSynthDots(s1);
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

if 0
tab1 = analyse(s, N, sprintf('%d_sdots.tif', s.nIter));
df_writeTable(tab1, 'tab1.csv');

tab2 = analyse(s, N, sprintf('dl2_%d_sdots.tif', s.nIter));
df_writeTable(tab2, 'tab2.csv');
end

if 1
s.nIter = 100;
tab3 = analyse(s, N, sprintf('dl2_%d_sdots.tif', s.nIter));
df_writeTable(tab3, 'tab3.csv');
end
if 0
s.nIter = 200;
tab4 = analyse(s, N, sprintf('dl2_%d_sdots.tif', s.nIter));
df_writeTable(tab4, 'tab4.csv');
end
if 0
s.nIter = 0;
tab5 = analyse(s, N, sprintf('sdots.tif', s.nIter));
df_writeTable(tab5, 'tab5.csv');
end
 
keyboard

% writetable(tab, 'tab.csv')
% save as table circumvent whimsical rounding



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