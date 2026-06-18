function DSCsample = DSCtools_setID(DSCsample, id_or_source)
% function DSCsample = DSCtools_setID(DSCsample, id)
% function DSCsample = DSCtools_addCP(DSCsample, 'by_source')
%
% Sets the id(s) to specified value(s) or selects a source for auto-generation
% 
% INPUT:   DSCsample --> DSCtools data structure of sample(s) as returned by DSCtools_readFile(s) - (e.g. pcm)
%                 id --> String or cell array of strings to be used as IDs. 
%        'by_source' --> if specified, an ID is created based on the selected source:
%                        'by_filename'    --> use filename as ID
%                        'by_hash'        --> use hash of filename as ID
%                        'by_description' --> generate ID from experimental description (sample, mass, rate)
%                                             if there are multiple with the same description, a unique number is added
%                        'by_enumeration' --> number in order of occurrence
%
% OUTPUT:    DSCsample --> DSCsample structure with updated ID fields
%
%
% Author:  Andreas Sommer  -- Jun2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%


% create IDs if not specified as strings
if istext(id_or_source) && startsWith(id_or_source, 'by_')
   switch lower(id_or_source)
      case {'by_enum', 'by_enumeration', 'by_number'}
         id = make_id_from_enumeration(DSCsample);
      case {'by_file', 'by_filename', 'by_filespec'}
         id = make_id_from_filename(DSCsample);
      case {'by_hash', 'by_filenamehash'}
         id = make_id_from_filenamehash(DSCsample);
      case {'by_desc', 'by_description'}
         id = make_id_from_description(DSCsample);
      otherwise
         error('Unknown source specified: %s', source);
   end
else
   id = id_or_source;
end

% from here on, we only use "id"
clear id_or_source

% if we have a single string, wrap it in a cell 
if istext(id)
   id = {id};
end
   
% ensure ID cell string is same size as number of DSC samples
nDSCsamples = length(DSCsample);
nIDs        = length(id);
if nDSCsamples ~= nIDs
   error('Number of IDs (%d) does not match number of DSC samples (%d)', nIDs, nDSCsamples);
end

% transfer IDs from cell string
for k = 1:nIDs
   DSCsample(k).ID = id{k};
end

% finito 
return

end


%% HELPERS

function id = make_id_from_enumeration(DSCsample)
   % just enumerate
   count = length(DSCsample);
   id = arrayfun(@num2str, 1:count, 'UniformOutput',false);
end

function id = make_id_from_filename(DSCsample)
   % use filename as ID
   count = length(DSCsample);
   id = cell(count, 1);
   for k = 1:count
      id{k} = DSCsample(k).rawData.fileSpec;
   end
end

function id = make_id_from_filenamehash(DSCsample)
   % use hash of filename as ID
   count = length(DSCsample);
   id = cell(count, 1);
   for k = 1:count
      id{k} = ADLER32(DSCsample(k).rawData.fileSpec);
   end
end

function id = make_id_from_description(DSCsample)
   % generate ID from sample, mass, rate
   count = length(DSCsample);
   id = cell(count, 1);
   for k = 1:count
      samp = DSCsample(k).rawData.desc.SAMPLE;
      mass = DSCsample(k).rawData.desc.SAMPLEMASSmg;
      rate = DSCsample(k).rawData.Tinfo.Tstep;
      id{k} = sprintf('%s__%smg__rate%+5.2f', samp, mass, rate);
   end
   % if there are multiple with the same ID now, add unique number to each
   [C, ~, ic] = unique(id, 'stable');
   if (length(C) ~= length(id))
      % occurrence number for each element (1st, 2nd, 3rd, ...)
      occVec = zeros(numel(strings), 1);
      for k = 1:numel(C)                  % run through unique strings
         idx = find(ic == k);             % all positions of this string
         occVec(idx) = 1:numel(idx);      % 1,2,3,... for this string
      end
      % append number to every ID
      for k = 1:length(id)
         id{k} = sprintf('%s__rep%d', id{k}, occVec(k));
      end
   end
end