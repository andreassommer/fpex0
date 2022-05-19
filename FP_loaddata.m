% FP_loaddata
% 
% Script that loads DSC measurement data 
%
%
% Copyright 2016-2022, Andreas Sommer  code@andreas-sommer.eu
%
% Copying and distribution of this file, with or without modification, are permitted in any medium without royalty, 
% provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.


% Global configuration
global FPEX0

% Baseline-Generation: enable/disable halt on warnings
global DSC204_getBaseline_haltOnWarning
DSC204_getBaseline_haltOnWarning = true;


% Select measurement files / file mask
withCorr = false;
if withCorr, corrmask = 'mit'; else corrmask = 'ohne'; end %#ok<*UNRCH,SEPEX>
measdir = FPEX0.paths.measurements;
dscdatafilemask   = fullfile(measdir, 'Messungen', sprintf('ExpDat_16-4*%sKorr*_H.csv', corrmask));
referencefilemask = fullfile(measdir, 'Waermekapazitaet_Saphirmessung', 'Sap-Kurve_10Kmin_H_Segment_7.csv');

% load dsc data
dscdata = DSC204_readFiles(dscdatafilemask);
dscdata = DSC204_sortByTstep(dscdata, 'ascend');

% DEBUG: remove 0,3K/min measurements (they contain unreliable temperatures with multiple same values)
idx03   = [dscdata.rate]==0.3 ;  % indices of 0.3 K/min measurements (testing for equality works here)
dscdata = dscdata(~idx03);       % exclude found indices

% DEBUG:  mitKorr: apply smoothing to ID 16-416, 16-417, 16-418 (corrected signals have high oscillations)
%        ohneKorr: apply smoothing to ID 16-416       (only that uncorrected signal has high oscillations)
if withCorr
   smoothIDs = {'16-416', '16-417', '16-418'};
else
   smoothIDs = {'16-416'};
end
for k=1:length(dscdata)
   if ismember(dscdata(k).desc.IDENTITY, smoothIDs) && (dscdata(k).rate == 0.6)
      fprintf('DEBUG: Applying smoothing to ID %s\n', dscdata(k).desc.IDENTITY);
      dscdata(k).data(:,3) = smoothdata(dscdata(k).data(:,3), 'sgolay', 75);
   end
end


% load saphire (reference) dsc data
dscsaphire = DSC204_readFiles(referencefilemask);
dscsaphire = DSC204_sortByTstep(dscsaphire, 'ascend');

% add cp curves
dscdata = DSC204_addCP(dscdata, dscsaphire);
