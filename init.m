% INIT
% INIT() Sets paths for FPEX0
%
%
% Copyright 2016-2022, Andreas Sommer  code@andreas-sommer.eu
%
% Copying and distribution of this file, with or without modification, are permitted in any medium without royalty, 
% provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.


% global base directory
global fpex0basedir
fpex0basedir = pwd();


% For working on a network share: (only on Windows)
if ispc
   disp('Windows detected. Invoking special policies for network shares.');
   if verLessThan('matlab', '9.3')
      system_dependent('RemoteCWDPolicy', 'Reload');
   end
   system_dependent('RemotePathPolicy', 'Reload');
end


% prepare paths
clear FPEX0paths
FPEX0paths.DSC204       = fullfile(fpex0basedir,'./DSC204');       % tools for DSC204 data processing
FPEX0paths.optionlists  = fullfile(fpex0basedir,'./optionlists');  % key-value pairs as option lists 
FPEX0paths.mmtools      = fullfile(fpex0basedir,'./mmtools');      % miscellaneous matlab tools
FPEX0paths.measurements = ''; % set here the paths to measurements, if necessary
% FPEX0paths.whatever   = fullfile(fpex0basedir,'./whatever');     % add whatever needed


% add the paths
for xTMPfieldname = fieldnames(FPEX0paths)'
   % extract full pathnames from the fields
   xTMPfieldname = xTMPfieldname{1};   %#ok<FXSET>  % xTMPfiledname is a 1-element cell array
   xTMPfullpath = FPEX0paths.(xTMPfieldname);   
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
   else % path does not exist: warn or error
      if isempty(xTMPfullpath)
         warning('Path for %s is missing. You should update init.m accordingly.', xTMPfieldname);
      else
         error('Missing path: %s',xTMPfullpath);
      end
   end
end


% cleanup
clear FPEX0paths
clear xTMPfullpath 
clear xTMPfieldname



 
