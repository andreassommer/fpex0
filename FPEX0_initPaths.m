% FPEX0_initPaths
% FPEX0_initPaths() Sets paths for FPEX0
%
% Adjust to your paths in the PATH PREPARATIONS section
%
% Copyright 2016-2022
% Andreas Sommer
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%
% Copying and distribution of this file, with or without modification, are permitted in any medium without royalty, 
% provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.


% where am I ?
[xTMPselfpath, ~, ~] = fileparts(which('FPEX0_initPaths.m'));

% access global path settings
global FPEX0paths

% setup
FPEX0paths.base = xTMPselfpath;  % global base directory


% For working on a network share: (only on Windows)
if ispc
   disp('Windows detected. Invoking special policies for network shares.');
   if verLessThan('matlab', '9.3')
      system_dependent('RemoteCWDPolicy', 'Reload');
   end
   system_dependent('RemotePathPolicy', 'Reload');
end



% PATHS PREPARATIONS
% ==================
FPEX0paths.DSC204       = fullfile(FPEX0paths.base, './DSC204');       % tools for DSC204 data processing
FPEX0paths.optionlists  = fullfile(FPEX0paths.base, './optionlists');  % key-value pairs as option lists 
FPEX0paths.mmtools      = fullfile(FPEX0paths.base, './mmtools');      % miscellaneous matlab tools
FPEX0paths.measurements = fullfile(FPEX0paths.base, '../DSC204_F1_Phoenix_Messungen'); % paths to measurements
%FPEX0paths.whatever     = fullfile(FPEX0config.base,'./whatever');     % add whatever path is ggggggggneeded



% transform the paths to absolute paths
for xTMPfieldname = fieldnames(FPEX0paths)'
   xTMPfullpath = what(FPEX0paths.(xTMPfieldname{1}));  % extract path info 
   FPEX0paths.(xTMPfieldname{1}) = xTMPfullpath.path;   % write back absolute path
end

% add the paths
for xTMPfieldname = fieldnames(FPEX0paths)'
   % extract full pathnames from the fields
   xTMPfieldname = xTMPfieldname{1};             %#ok<FXSET>  % xTMPfiledname is a 1-element cell array
   xTMPfullpath = FPEX0paths.(xTMPfieldname);   % extract path name (may be relative)
   % check if directory exists
   if exist(xTMPfullpath,'dir')
      fprintf('Processing %s \t',xTMPfullpath);
      % only add path if not yet in path
      if ~ismember(xTMPfullpath, strsplit(path, pathsep))
         fprintf('-- Adding to path.\n');
         addpath(xTMPfullpath);
      else
         fprintf('-- Already in path, skipping.\n');
      end
   else % path does not exist: warn
      if isempty(xTMPfullpath)
         warning('Path for %s is missing. You should update init.m accordingly.', xTMPfieldname);
      else
         warning('Path for "%s" does not exist: %s', xTMPfieldname, xTMPfullpath);
      end
   end
end


% activate the DSC204 settings
DSC204_settings('ACTIVATE');


% cleanup
clear xTMPfullpath 
clear xTMPfieldname
clear xTMPselfpath


 
