dirs = dir('*');

for kk = 1:numel(dirs)
    folder = dirs(kk).name;
    if isdir(folder) && numel(folder) > 2
        disp(folder)
        files = dir([folder filesep() '*tif']);
        v = imread([folder filesep() files(1).name]);
        V = zeros([size(v), numel(files)]);
        for ll = 1:numel(files)
            V(:,:,ll) = imread([folder filesep() files(ll).name]);
        end
        outname = sprintf('%s.tif', folder);
        disp(outname)        
        df_writeTif(uint16(V), outname);
    end
end