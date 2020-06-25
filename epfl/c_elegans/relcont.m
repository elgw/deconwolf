function relcont()
addpath('../');
files(1).name = 'CElegans-FITC.tif';
files(1).tag = 'Original';
files(end+1).name = 'dw_CElegans-FITC.tif';
files(end).tag = 'dw';
files(end+1).name = 'dw0_CElegans-FITC.tif';
files(end).tag = 'dw0';
files(end+1).name = 'dw1_CElegans-FITC.tif';
files(end).tag = 'dw1';
files(end+1).name = 'Final Display of RL_FITC.tif';
files(end).tag = 'deLa2';

figure,
hold on

legendstrs = {};
for kk = 1:numel(files)
    file = files(kk).name;
    I = double(df_readTif(file));
    
    S = I(:,:,48);
    p0 = [247, 255];
    p1 = [215, 413];
    len = norm(p0-p1);
    x = linspace(p0(1), p1(1), round(len));
    y = linspace(p0(2), p1(2), round(len));
    
    curve = interpn(S, x, y);
    cmax = max(curve);
    cmin = min(curve(40:120));
    rc = 1-cmin/cmax;
    files(kk).relContrast = rc;
    files(kk).mem_Mb = dwGetMem(file)/1000;
    files(kk).time_s = dwGetTime(file);
    plot(curve/max(curve))
        
    fprintf('File: %s, rel contrast: %f\n', file, rc);
    %volumeSlide(I)
    %keyboard        
    legendstrs{end+1} = files(kk).tag;
end

legend(legendstrs, 'interpreter', 'none')
xlabel('Position on the line [pixels]')
ylabel('Scaled intensity')
dprintpdf('relcont_c_elegans', 'publish', 'w', 20)

tab = struct2table(files);
fprintf('Writing results to tab.csv');
df_writeTable(tab, 'tab.csv');

end