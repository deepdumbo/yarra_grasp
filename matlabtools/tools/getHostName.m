% Function for determining what machine Matlab is currently running on
% 
% Written by Marnix Maas (Marnix.Maas@Radboudumc.nl), 1-oct-2018

function [hostname, os] = getHostName()
os = '';
hostname = '';
if ispc
    % We're on Windows
    os = 'Windows';
    [~, hostname] = system('hostname');
    hostname = hostname(1:end-1);
elseif isunix
    % We're on Unix/Linux
    os = 'Linux';
    [~, hostname] = system('uname -n');
    hostname = hostname(1:end-1);
else
    % Not implemented for anything else than Win or Ux
    warning('Hostname detection not implemented for operating systems other than Windows or Unix/Linux');
end
