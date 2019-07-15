function [twix, rawdata] = loadData(fileName, varargin)
% Load Twix header data and k-space data from raw data file specified by
% filename. 
% Optional arguments can be specified:
% 
% Written by Marnix Maas, RadboudUMC (Marnix.Maas@radboudumc.nl)
% September 2018

% Parse input arguments
% p = parseInputs(fileName, varargin);
p = parseInputs(varargin);
parse(p, fileName, varargin{:});

doLoadTwixFile  = p.Results.doLoadTwixFile;
doSaveTwixFile  = p.Results.doSaveTwixFile;
doUnsorted      = p.Results.doUnsorted;
twixFilePath    = p.Results.twixFilePath;
twixFileName    = p.Results.twixFileName;
removeOS        = p.Results.removeOS;
ignoreSeg       = p.Results.ignoreSeg;
channels        = p.Results.channels;
spokes          = p.Results.spokes;
partitions      = p.Results.partitions;
echoes          = p.Results.echoes;
repeats         = p.Results.repeats;

% ## Read metadata from the Siemens TWIX file
fullPath = fullfile(twixFilePath, twixFileName);
if doLoadTwixFile
    % Load header info from previously saved file: faster
%     if isfile(fullPath)   % only works from R2017b
    if exist(fullPath, 'file')
		load(fullPath)
    else
        msg = 'No valid file name for Twix data provided: aborting';
%         logRecon(msg, fullfile(temp_path,pars.logFileName, pars.doShowLogMsg));
        error(msg)
    end
else
    twix = mapVBVD(fileName);
    
    % For VD software, the TWIX file may contain adjustment data. We only want
    % to look at the image data for now.
    % Depending on how the data was transferred (e.g. via Twix file copy or
    % Yarra), the twix data object may be a struct or a cell array.
    if isa(twix, 'cell')
        twix = twix{2};     % MCM: should this be 'end' instead of '2'?
    end
    
    % Save twix file for quicker loading later
    if doSaveTwixFile
        save(strcat(fullPath), 'twix')
    end
end

% Set some parameters for reading the raw data
twix.image.flagRemoveOS     = removeOS;     % Remove readout oversampling to reduce memory load
twix.image.flagIgnoreSeg    = ignoreSeg;    % Merge k-space 'segments' together to reduce memory load

if channels == 0
    channels = 1:twix.image.NCha;
end
if spokes == 0
    spokes = 1:twix.image.NLin;
%     spokes = 1:twix.image.NLinMeas;         % .NLin doesn't exist in 3d radial dataset... check if this one works for 3D SoS data
end
if partitions == 0
    partitions = 1:twix.image.NPar;
end
if echoes == 0
    echoes = 1:twix.image.NEco;
end
if repeats == 0
    repeats = 1:twix.image.NRep;
end

% Get the k-space data. 
if doUnsorted
    rawdata = twix.image.readTW(length(spokes)*length(repeats));
    rawdata = squeeze(rawdata(:,channels,:));
else
    % Data comes in as [samples,channels,spokes,partitions,~,~,~,echoes,repeats]
    rawdata = squeeze(twix.image(:,channels,spokes,partitions,:,:,:,echoes,repeats));
end

end

 


% function p = parseInputs(fileName, varargin)
function p = parseInputs(varargin)
p = inputParser;
% default values
defDoLoadTwixFile   = 0;
defDoSaveTwixFile   = 1;
defDoUnsorted       = 0;
defTwixFileName     = '';
defTwixFilePath 	= '../../../twixdata';
defRemoveOS         = 0;
defIgnoreSeg        = 1;
defChannels         = 0;
defSpokes           = 0;
defPartitions       = 0;
defEchoes           = 0;
defRepeats          = 0;

% fill parser object
addRequired(p,'fileName',@ischar);
addParameter(p,'doLoadTwixFile',defDoLoadTwixFile,@isnumeric);
addParameter(p,'doSaveTwixFile',defDoSaveTwixFile,@isnumeric);
addParameter(p,'doUnsorted',defDoUnsorted,@isnumeric);
addParameter(p,'twixFilePath',defTwixFilePath,@ischar);
addParameter(p,'twixFileName',defTwixFileName,@ischar);
addParameter(p,'removeOS',defRemoveOS,@isnumeric);
addParameter(p,'ignoreSeg',defIgnoreSeg,@isnumeric);
addParameter(p,'channels',defChannels,@isnumeric);
addParameter(p,'spokes',defSpokes,@isnumeric);
addParameter(p,'partitions',defPartitions,@isnumeric);
addParameter(p,'echoes',defEchoes,@isnumeric);
addParameter(p,'repeats',defRepeats,@isnumeric);
end