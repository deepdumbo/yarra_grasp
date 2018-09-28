function [twix, rawdata] = loadDataGoldenAngle(fileName, varargin)
% Load Twix header data and k-space data from raw data file specified by
% filename. 
% Optional arguments can be specified:
% 
% Written by Marnix Maas, RadboudUMC (Marnix.Maas@radboudumc.nl)
% September 2018

% Parse input arguments
p = parseInputs(fileName, varargin);
parse(p, fileName, varargin{:});

doLoadTwixFile  = p.Results.doLoadTwixFile;
doSaveTwixFile  = p.Results.doSaveTwixFile;
twixFilePath    = p.Results.twixFilePath;
twixFileName    = p.Results.twixFileName;
removeOS        = p.Results.removeOS;
channels        = p.Results.channels;
spokes          = p.Results.spokes;
partitions      = p.Results.partitions;
echoes          = p.Results.echoes;

% ## Read metadata from the Siemens TWIX file
fullPath = fullfile(twixFilePath, twixFileName);
if doLoadTwixFile
    % Load header info from previously saved file: faster
%     if isfile(fullPath)   % only works from R2017b
    if exist(fullPath, 'file')
		load(fullPath)
    else
        msg = 'No valid file name for Twix data provided: aborting';
        error(msg)
    end
else
    twix = mapVBVD(fileName);
    
    % For VD software, the TWIX file may contain adjustment data. We only want
    % to look at the image data for now.
    twix = twix{2};
    if doSaveTwixFile
        save(strcat(fullPath), 'twix')
    end
end

% Set some parameters for reading the raw data
twix.image.flagRemoveOS = removeOS;         % Remove readout oversampling to reduce memory load

if channels == 0
    channels = 1:twix.image.NCha;
end
if spokes == 0
    spokes = 1:twix.image.NLin;
end
if partitions == 0
    partitions = 1:twix.image.NPar;
end
if echoes == 0
    echoes = 1:twix.image.NEco;
end

% Get the k-space data. Data comes in as [samples,channels,spokes,partitions,~,~,~,echoes]
rawdata = squeeze(twix.image(:,channels,spokes,partitions,:,:,:,echoes));       % add selection of echoes here

end

 


function p = parseInputs(fileName, varargin)
p = inputParser;
% default values
defDoLoadTwixFile   = 0;
defDoSaveTwixFile   = 1;
defTwixFileName     = '';
defTwixFilePath 	= '../../../twixdata';
defRemoveOS         = 0;
defChannels         = 0;
defSpokes           = 0;
defPartitions       = 0;
defEchoes           = 0;

% fill parser object
addRequired(p,'fileName',@ischar);
addParameter(p,'doLoadTwixFile',defDoLoadTwixFile,@isnumeric);
addParameter(p,'doSaveTwixFile',defDoSaveTwixFile,@isnumeric);
addParameter(p,'twixFilePath',defTwixFilePath,@ischar);
addParameter(p,'twixFileName',defTwixFileName,@ischar);
addParameter(p,'removeOS',defRemoveOS,@isnumeric);
addParameter(p,'channels',defChannels,@isnumeric);
addParameter(p,'spokes',defSpokes,@isnumeric);
addParameter(p,'partitions',defPartitions,@isnumeric);
addParameter(p,'echoes',defEchoes,@isnumeric);
end