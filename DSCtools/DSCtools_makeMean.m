function DSCmean = DSCtools_makeMean(DSCdata, Trange)
% DSCmean = DSCtools_makeMean(DSCdata, Trange)
%
% Calculates the mean value of the measurements in DSCdata and returns a "mean measurement" dsc.
% This is useful for analysis, rather not for further processing.
% 
% INPUT:  DSCdata --> DSC data structure as returned by DSCtools_readFile()
%          Trange --> [Tmin Tmax] temperature range to be used
%
% OUTPUT: DSCmean --> artificial DSC data structure with mean values from DSCdata.
%
% NOTE:  Some of the internal data are set to artificial texts (like ID or fileSpec), 
%        others are removed (like fileHeader or rawData.datMat) .
%
% Author: Andreas Sommer, Jun2026
% email@andreas-sommer.eu


% get prototype and remove things that we better do not touch
DSCmean = DSCdata(1);
DSCmean.rawData.fileHeader = {};
DSCmean.rawData.fileSpec = '--MEAN--';
DSCmean.rawData.datMat = [];
DSCmean.ID = '--MEAN--';

% it does not make sense, but if the user puts in data from mixed heat rates, we also manipulate it
rates = [DSCdata.rate];    DSCmean.rate = mean(rates);
masses = [DSCdata.mass];   DSCmean.mass = mean(masses);

% retrieve common range of the temperatures
Tmin =  max( arrayfun(@(x)    min(x.data.T), DSCdata) );
Tmax =  min( arrayfun(@(x)    max(x.data.T), DSCdata) );
Tn   = mean( arrayfun(@(x) length(x.data.T), DSCdata) );

% check user requested temperature range
if Trange(1) < Tmin, warning('Minimum Temperature must be >= %g', Tmin), Trange(1) = Tmin; end
if Trange(2) > Tmax, warning('Maximum Temperature must be <= %g', Tmax), Trange(2) = Tmax; end
Tmin = Trange(1);
Tmax = Trange(2);
T = linspace(Tmin, Tmax, Tn);

% interpolate the data fields and take the mean
DSCmean.data.T  = T;
DSCmean.data.uV = calcMean('uV', T,  DSCdata );
DSCmean.data.mW = calcMean('mW', T,  DSCdata );
DSCmean.data.sf = calcMean('sf', T,  DSCdata );
DSCmean.data.t  = [0  (T / DSCmean.rate) ];

% possibly the CP field
if any(arrayfun(@(x) isfield(x, 'cp'), DSCdata))
   warning('CP mean computation not yet implemented')
end

% setup Tinfo field
DSCmean.Tinfo.Tstart = T(1);
DSCmean.Tinfo.Tend   = T(end);
DSCmean.Tinfo.Tstep  = DSCdata.rate;

% finito
return


end % of function

%% HELPERS
function Yq = interpolateData(T,Y,Tq)
   idx = (T>=Tq(1) & T<=Tq(end));
   Yq = interp1(T(idx), Y(idx), Tq, 'linear');
end

function Y = calcMean(field, T, dsc)
  Y = arrayfun( @(x) interpolateData(x.data.T, x.data.(field), T), dsc, 'UniformOutput', false );  % get data in cell array
  Y = cell2mat( reshape(Y, [], 1) );        % put cells one beneath the other
  Y = mean(Y);                              % compute column-wise mean
end