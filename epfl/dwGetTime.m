function t = dwGetTime(file)

logFile = [file, '.log.txt'];
fid = fopen(logFile, "r");
if fid == -1
    t = 0;
    return;
end
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if strfind(tline, 'Took:') == 1
        t = str2num(tline(6:end-2));
    end
end
fclose(fid);
end
