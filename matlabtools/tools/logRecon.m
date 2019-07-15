function logRecon(msg, fileName, show)
% Write messages generated during reconstruction to a log file. Optionally
% display them on screen as well.
% Inputs:
% - msg:        message to be written/displayed
% - fileName:   full path to log file
% - disp:       flag for displaying messages on screen

fid = fopen(fileName, 'a');
if fid == -1
  error('Cannot open log file.');
end
% hostname = getHostName();
[~, hostname] = system('hostname');
st = dbstack;
caller = st(2).name;
fprintf(fid, '%s: %s: %s: %s\n', datestr(now, 0), hostname, caller, msg);
fclose(fid);
if show
    disp(msg)
end