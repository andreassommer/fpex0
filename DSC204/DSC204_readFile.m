function dataStruct = DSC204_readFile(fileSpec)
   % function data = DSC204_readFile(fileSpec)
   %
   % Reads DSC204 CSV file and stores everything in data structure.
   %
   % INPUT:  fileSpec --> string describing the file path and name
   %
   % OUTPUT:  structure with fields:
   %           .data: [8160x4 double]                                              % read data
   %  .columnHeaders: {'Temp./'C'  'Time/min'  'DSC/(uV/mg)'  'Sensit./(uV/mW)'}   % column headers, where the o in oC marks degree 
   %     .fileHeader: {1x32 cell}                                                  % file headers as strings
   %       .fileSpec: '~/DSCMessungen/407/ExpDat_16-407-3_ohneKorr_0,3Kmin_H.csv'  % file spec of read file
   %        .fileEnc: [1x1 struct]                                                 % file encoding
   %           .desc: [1x1 struct]                                                 % description information (read from file headers)
   %          .Tinfo: [1x1 struct]                                                 % temperature profile (read from file headers)
   %           .rate: 1.25                                                         % heating/cooling rate (same as in .Tinfo.Tstep)
   %           .mass: 10.47                                                        % sample mass
   %  .multiplescans: 0                                                            % flag indicating if multiple data are read
   %
   %
   % NOTE: I don't have any documentation about the file format. Everything is done by reverse-engineering.
   %       Use at your own risk and proceed with fingers crossed.
   %
   % Author:  Andreas Sommer, Mar2017, Aug2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %

 
   % default encodings
   enc = DSC204_getEncodingDefaults();

   % open file
   fid = fopen(fileSpec, 'r', enc.machineformat, enc.encoding);
 
   % read the header (possibly updates the encoding)
   [desc, enc, fileHeaderStrings] = DSC204_readFileHeader(fid, enc);
   
   % try to extract IDENTITY from dsc, otherwise create one from filename
   if isfield(desc,'IDENTITY')
      IDENTITY = desc.IDENTITY;
   else
      IDENTITY = ADLER32(fileSpec);
   end
   
   % retrieve the column headers
   columnHeaders = DSC204_readColumnHeaders(fid, enc);
   
   % read data
   dataMat = DSC204_readData(fid, enc, length(columnHeaders), NaN);
   
   % close file
   fclose(fid);
   
   % extract information from description
   Tinfo = extractTinfo(desc);       % temperature driving information
   mass  = extractSampleMass(desc);  % sample mass
   if ~isempty(Tinfo)
      rate = Tinfo.Tstep;
   else
      rate = []; 
   end  % heating/cooling rate
   
   % check if there were multiple DSC scans in one file
   multiplescans = (length(desc) > 1);
      
   % store the raw data in separate substructure
   rawData = struct(...
                        'columnHeaders', {columnHeaders}       ,...
                        'fileHeader'   , {fileHeaderStrings}   ,...
                        'fileSpec'     , fileSpec              ,...
                        'fileEnc'      , enc                   ,...
                        'dataMat'      , dataMat               ,...
                        'desc'         , desc                  ,...
                        'Tinfo'        , Tinfo                 ,...
                        'multiplescans', multiplescans          ...
                        );

   % generate result structure
   dataStruct  = struct(...
                        'ID'           , IDENTITY              ,...
                        'rawData'      , rawData               ,...
                        'data'         , []                    ,...   % gets filled later
                        'rate'         , rate                  ,...
                        'mass'         , mass                   ...
                        );

   % add quick accessors
   dataStruct = DSC204_addQuickAccessors(dataStruct, -inf, +inf);
   
   % finito
   return
   
   
   % ===========================================
   % helpers
    

   function Tinfo = extractTinfo(desc)
      try
         for k = 1:length(desc)
            Tinfo(k) = DSC204_extractRANGE(desc(k).RANGE, enc);                 %#ok<AGROW>
         end
      catch
         warning('Unable to find temperature driving information.')
         Tinfo = [];
      end
   end     

   function mass = extractSampleMass(desc)
      mass = zeros(1, length(desc));
      for k = 1:length(desc)
         try
            massString = desc(k).SAMPLEMASSmg;
            massString = DSC204_replaceDecimalDelimiter(massString, enc);
            [m, success] = str2num(massString);  % ATTENTION: str2num understands '323.3d0', str2double does not.
            if ~success, error('Can''t convert %s to double.', massString), end
            mass(k) = m;
         catch err
            warning('Error while searching sample mass: %s.', err.message)
         end
      end
      
   end
   
                    
end

   