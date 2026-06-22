function dscdata = DSCtools_readFiles(fileSpecs)
   % function dscdata = DSCtools_readFiles(fileSpecs)
   %
   % Reads specified DSCtools CSV files and stores everything in dataStruct.
   %
   % INPUT:     fileSpec --> 1) cell array of strings describing the file paths and names
   %                         2) a single string describing the files to load in wildcard form (e.g. 'DSC*.csv');
   %
   % OUTPUT:     dscdata --> Struct array with contents as read by DSCtools_readFile.
   %                         See there for documentation.
   %
   %
   % Author: Andreas Sommer  --  Mar2017, Aug2022, Jun2026
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   
   % if fileSpecs is a string, read files
   if ischar(fileSpecs)
      [pathSpec, ~, ~] = fileparts(fileSpecs);
      dirlist  = dir(fileSpecs);
      if isempty(dirlist)
         error('Could not find any file matching "%s"', fileSpecs)
      end
      fileSpecs = fullfile(pathSpec,{dirlist.name});
   end
   
   % walk through the file specs
   count = length(fileSpecs);
   dataStructs = cell(count,1);
   for k = 1:count
      fileSpec = fileSpecs{k};
      try
         makeMessage('Reading %s  [%d/%d]\n', fileSpec, k, count);
         contents = DSCtools_readFile(fileSpec);
      catch err
         warning('\nWhile reading %s, caught error: %s. SKIPPING FILE!', fileSpec, err.message);
         contents = [];
      end
      dataStructs{k} = contents;
   end
   
   % transfer them into single array
   nonEmptyIdx = ~cellfun(@isempty, dataStructs);
   dscdata  = cell2mat(dataStructs(nonEmptyIdx));
   
   
end