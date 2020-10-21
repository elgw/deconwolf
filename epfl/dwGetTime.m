function t = dwGetTime(file)
t = -1;

logFile = [file, '.log.txt'];
fid = fopen(logFile, "r");
if fid == -1   
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
