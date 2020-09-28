function t = dwGetMem(file)
logFile = [file, '.log.txt'];
fid = fopen(logFile, "r");
if fid == -1
    t = 0;
    return;
end
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if contains(tline, 'peakMemory:') == 1
        tline = strsplit(tline, ' ');
        t = str2num(tline{2});
    end
end
fclose(fid);
end

