function varargout = DSC204test_makeMeasurementsMatFile(datafilemask, referencefilemask, outputfile)
   % measobj = DSC204test_makeMeasurementsMatFile(datafilemask, referencefilemask, outputfile)
   %
   % Creates a measurement file that contains the dsc measurement data as well as a measurement function of
   % two variables (heating rate and temperature).
   %
   % INPUT:      datafilemask --> DSC measurement files to use
   %        referencefilemask --> DSC measurement reference files to use (sapphire)
   %               outputfile --> file to write to
   %
   % OUTPUT:    measobj --> object as written to file
   %
   %
   % Usage example:
   %
   % mdata = load('DSC_FP_measurements.mat');  % loads the data (the file contains a single structure called 'mdata')
   % mdata = mdata.mdata;                      % reduce one level of nesting
   % mtime  = 0.6;                             % time point of queried measurements
   % mspace = 30:0.1:220;                      % (vector of) spacial coordinates of queries measurements
   % ttol = 0.2;                               % tolerance of finding mtime (default: 0.1)
   %
   % X = mdata.mfun(mtime, mspace, ttol);      % retrieve data
   % plot(mspace,X,'.-')                       % and plot it
   %
   %
   % Structure of mdata:
   %
   % mdata.tvec --> vector of measurment times    (i.e. heating rates)
   % mdata.mfun --> handle to measurement function
   % 
   
   
   % default file mask and output file
   defaultpath = '/home/asommer/Documents/modELTES/DSC204_F1_Phoenix_Messungen/alles/';
   if (nargin < 1) || isempty(datafilemask)
      datafilemask = fullfile(defaultpath, 'ExpDat_16-408*mitKorr*_H.csv');
   end
   if (nargin < 2) || isempty(referencefilemask)
      referencefilemask = fullfile(defaultpath, 'Sap-Kurve_10Kmin_H_Segment_7.csv');
   end
   if (nargin < 3)
      outputfile = 'DSC_FP_measurements.mat'; 
   end

   % load dsc data and sort by rate
   dscdata = DSC204_readFiles(datafilemask);
   dscdata = DSC204_sortByTstep(dscdata, 'ascend');
   
   % load saphire (reference) dsc data
   dscsaphire = DSC204_readFiles(referencefilemask);
   
   % add cp curves
   dsc = DSC204_addCP(dscdata, dscsaphire);
   
   % heating rates = measurement time points in PF
   rates = [dsc.rate]; 
   
   % extract cp values and make interpolation function
   for k = 1:length(dsc)
      T{k}       = dsc(k).cp.T;           %#ok<*AGROW>
      lcpdata{k} = dsc(k).cp.latentdata;  % lcp = latent cp   
   end
   
   % make object
   mdata.tvec     = rates;   % measurement times in FP  =  heating rates
   mdata.T        = T;       % raw data time points
   mdata.latentcp = lcpdata; % raw data 
   
   % store the dsc object in outputfile
   if ~isempty(outputfile)
      save(outputfile, 'mdata', 'datafilemask', 'referencefilemask');
      fprintf('\n\nWritten file: %s\n', outputfile);
   end

   % set output variable
   if (nargout > 0)
      varargout{1} = mdata;
   end
   
end