function alldata = DSC204_loadall(filemask)
  % Loads all DSC measurement files in specified filemask (in alphabetical order)
  %
  % INPUT:  filemask --> filemask for files to be loaded, e.g. '*ohneKorr*_H.csv'
  %
  % OUTPUT:  alldata --> structure array containing the results from DSC204_readFile
  %                      applied to each file
  %
  %
  % Author: Andreas Sommer, Mar2017
  % andreas.sommer@iwr.uni-heidelberg.de
  % email@andreas-sommer.eu
  %
  
  warning('DSC204_loadall is deprecated. Use DSC204_readFiles instead.')
  
  % no dir specified? use default
  if nargin<1
     filemask = fullfile('/home/asommer/Projekte/modELTES/DSC204_F1_Phoenix_Messungen/Messungen/alle','*.csv');
  end
  
  % extract path from filemask
  [pathstr, ~, ~] = fileparts(filemask); 
  
  % get csv file list, sort it alphabetically
  files = dir(fullfile(filemask));
  [~, sortidx] = sort({files(:).name});
  files = files(sortidx);
  
  % walk through files
  for k = 1:length(files)
     % get current file
     fstruct = files(k);
     % skip directory (who names a directory with ending .csv?)
     if fstruct.isdir, continue, end;
     % load the file
     fprintf('Loading: %s  \t', fstruct.name); 
     tic
     alldata(k) = DSC204_readFile(fullfile(pathstr,fstruct.name));
     toc
  end
  fprintf('\n\n');
  
  
end