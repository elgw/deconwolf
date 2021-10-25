function getmse()

psf = df_readTif('PSF-Bars-stack.tif');

files(1).Name = 'Bars-G10-P15-stack.tif';
files(end).Tag = 'Synthetic SNR15';
files(end+1).Name = 'Bars-G10-P30-stack.tif';
files(end).Tag = 'Synthetic SNR30';

files(end+1).Name = 'dw50_Bars-G10-P15-stack.tif';
files(end).Tag = 'dw SNR15';
files(end+1).Name = 'dw50_Bars-G10-P30-stack.tif';
files(end).Tag = 'dw SNR30';

files(end+1).Name = 'dw50.bq1_Bars-G10-P15-stack.tif';
files(end).Tag = 'dw50.bq1 SNR15';
files(end+1).Name = 'dw50.bq1_Bars-G10-P30-stack.tif';
files(end).Tag = 'dw50.bq1 SNR30';

files(end+1).Name = 'dw50.bq0_Bars-G10-P15-stack.tif';
files(end).Tag = 'dw50.bq0 SNR15';
files(end+1).Name = 'dw50.bq0_Bars-G10-P30-stack.tif';
files(end).Tag = 'dw50.bq0 SNR30';

files(end+1).Name = 'P15_deconvolution_lab50.tif';
files(end).Tag = 'DL2RL50 SNR15';
files(end+1).Name = 'P30_deconvolution_lab50.tif';
files(end).Tag = 'DL2RL50 SNR30';

% We don't scale the reference image
Iref = df_readTif('Bars-stack.tif');
Iref = double(Iref); 

for ff = 1:numel(files)
    file = files(ff).Name;
    tag = files(ff).Tag;
    I = df_readTif(file);
    I = double(I);    
    I = I*mean(Iref(:))/mean(I(:)); % For easier opimization
    [scaling, mse] = fminunc(@(x) mseval(Iref, I*x), 0.1);
    fprintf('mse(%s,ref) = %f\n', tag, mse);
    files(ff).mse = mse;
    files(ff).time_s = dwGetTime(file);
    files(ff).mem_kb = dwGetMem(file);
end

tab = struct2table(files);
% writetable(tab, 'tab.csv')
% save as table circumvent whimsical rounding
df_writeTable(tab, 'bars50tab.csv');
end


function v = mseval(A, B)
 v = sqrt(mean((A(:)-B(:)).^2));
end

