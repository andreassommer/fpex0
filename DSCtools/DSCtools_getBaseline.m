function [blobj, bldata] = DSCtools_getBaseline(T, Y, bltype, blds_preset)
% [blobj, bldata] = DSCtools_getBaseline(T, Y, bltype, blds_preset)
%
% Retrieves various baseline function for specified DSCTOOLS_GLOBAL_SETTINGS.
%
%  INPUT:   T --> vector of x-values (e.g. temperatures)
%           Y --> vector of y-values (e.g. cp-values)
%      bltype --> type of baseline:  'linrange', 'ToTe', 'simple'
%                 add the suffix "_sig" to get the sigmoidal version   [default: '' to take from blds]
%        blds --> string describing the BaseLine Detection Setting     [default: '' to take from blds]
%
%  OUTPUT:  blobj --> baseline object of type DSCtools_Baseline, containing baseline function, onset, endset
%          bldata --> structure containing information about various baselines
%
%  NOTE: The following 6 baseline types are computed and stored in bldata
%          bl_linrange --> baseline built from detected linear segments
%          bl_ToTe     --> baseline built from linear segments defined by Tonset and Toffset
%          bl_simple   --> baseline built from detected linear segments
%          bl_*_sig    --> sigmoidal versions of the above
%
% Author:  Andreas Sommer, Apr2017, Jul2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%

% halt on warning?
global DSCTOOLS_GLOBAL_SETTINGS  %#ok<GVMIS>
haltOnWarning = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'getBaseline_haltOnWarning');

% retrieve load baseline detection settings from default if not specified
if (nargin < 4) || isempty(blds_preset)
   blds_preset = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'blds_preset');
end

% load baseline detection settings 
blds = DSCtools_getBaselineDetectionSettings(blds_preset);

% specify baseline type
if (nargin < 3) || isempty(bltype)
   bltype = blds.baselinetype;
end

% detect linear ranges according to baseline detection settings
makeMessage('Detecting linear ranges using preset "%s" ... ', blds_preset);
[~, peakPos] = max(Y); 
initlen_L       = floor(blds.L.initfraction * peakPos);
initlen_R       = floor(blds.R.initfraction * (length(T)-peakPos));
detectionArgs_L = {blds.L.reldevA, blds.L.reldevB, blds.L.reldevS2, blds.L.absdevA, blds.L.absdevB, blds.L.absdevS2};
detectionArgs_R = {blds.R.reldevA, blds.R.reldevB, blds.R.reldevS2, blds.R.absdevA, blds.R.absdevB, blds.R.absdevS2};
[idx_L, reg_L]  = DSCtools_detectLinearRange(T,Y,'left' ,initlen_L,detectionArgs_L{:});
[idx_R, reg_R]  = DSCtools_detectLinearRange(T,Y,'right',initlen_R,detectionArgs_R{:});
makeMessage('#raw', 'done.\n')

% make baselevel functions from linear ranges
bl_linrange     = DSCtools_getBaselinePrimitive(T, Y, idx_L, reg_L.a, reg_L.b, idx_R, reg_R.a, reg_R.b, 'linear'); 
bl_linrange_sig = DSCtools_getBaselinePrimitive(T, Y, idx_L, reg_L.a, reg_L.b, idx_R, reg_R.a, reg_R.b, 'sigmoidal'); 

% Estimate onset/endset: Point X where the data (X,Y) first falls below the LINEAR baseline (seen from peak maximum)
% with fallback using the ind_L or idx_R from linear parts
bloffset = blds.bloffset;
idxOnset  = find( Y(  1:peakPos) - bl_linrange.eval(T(  1:peakPos)) < bloffset, 1, 'last' );
idxEndset = find( Y(peakPos:end) - bl_linrange.eval(T(peakPos:end)) < bloffset, 1, 'first') + peakPos - 1;
if isempty(idxOnset),  idxOnset  = idx_L; end
if isempty(idxEndset), idxEndset = idx_R; end
Tonset  = T(idxOnset );
Tendset = T(idxEndset);

% Determine means of first and of last cp-values over a temperature window of 1 degree Celsius
Twindow = 1.0;  % window length w: one degree celsius
w = find(T-T(1)>=Twindow, 1, 'first');
mean_cp_L = mean(Y(1:w));
mean_cp_R = mean(Y(end-w:end));

% Determine baseline from Tonset and Tendset -- assuming linear parts left and right over whole temperature range
[slopeL, iceptL] = calc_mx_plus_b(T( 1 ), mean_cp_L, Tonset , Y(idxOnset) );
[slopeR, iceptR] = calc_mx_plus_b(T(end), mean_cp_R, Tendset, Y(idxEndset));
bl_ToTe     = DSCtools_getBaselinePrimitive(T, Y, idxOnset, iceptL, slopeL, idxEndset, iceptR, slopeR, 'linear');
bl_ToTe_sig = DSCtools_getBaselinePrimitive(T, Y, idxOnset, iceptL, slopeL, idxEndset, iceptR, slopeR, 'sigmoidal');


% Determine most simple baseline over whole temperature range
[slopeS, iceptS] = calc_mx_plus_b(T( 1 ), mean_cp_L, T(end) , mean_cp_R);
bl_simple     = DSCtools_getBaselinePrimitive(T, Y, 1, iceptS, slopeS, length(T), iceptS, slopeS, 'linear');
bl_simple_sig = DSCtools_getBaselinePrimitive(T, Y, 1, iceptS, slopeS, length(T), iceptS, slopeS, 'sigmoidal');


% small plausibility check:  slope of linear parts should be small
maxslope = 0.1;
if (abs(reg_L.b) > maxslope) || (abs(reg_R.b) > maxslope)
   makeMessage('#raw','\n');
   warning('Slope of linear part is large: left: %g, right: %g. There''s probably something wrong. Using proposed baseline instead!', reg_L.b, reg_R.b)
   oldbl = bl_linrange;
   if (abs(reg_L.b) > maxslope); reg_L.b = 0; reg_L.a = Y(1)  ; idx_L = 2          ; end
   if (abs(reg_R.b) > maxslope); reg_R.b = 0; reg_R.a = Y(end); idx_R = length(Y)-1; end
   bl_linrange = DSCtools_getBaselinePrimitive(T, Y, idx_L, reg_L.a, reg_L.b, idx_R, reg_R.a, reg_R.b, bltype);
   if haltOnWarning
      figure(271); clf; plot(T,Y,T,oldbl.eval(T),T,bl_linrange.eval(T)); 
      legend('data','detected baseline','proposed baseline','Location','best');
      makeMessage('Entering debugger... Press F5 to proceed.\n');
      keyboard  %#ok<KEYBOARDFUN>
   end
end

% DEBUG MODE:
%figure(271); clf; plot(X,Y,X,blfun(X)); legend('data', 'using baseline', 'Location','best'); disp('Entering debugger... Press F5 to proceed.'); keyboard

% collect baseline functions - first linear then sigmoidal
baselines.bl_linrange     = bl_linrange;
baselines.bl_ToTe         = bl_ToTe;
baselines.bl_simple       = bl_simple;
baselines.bl_linrange_sig = bl_linrange_sig;
baselines.bl_ToTe_sig     = bl_ToTe_sig;
baselines.bl_simple_sig   = bl_simple_sig;

% prepare output struct
bldata = struct();
bldata.reg_L  = reg_L;       % regression information left
bldata.reg_R  = reg_R;       % regression information right
bldata.linL   = T(idx_L);    % classic onset/endset estimation as ..
bldata.linR   = T(idx_R);    % .. end and start of linear parts (poor estimation)
bldata.onset  = Tonset;      % estimate of onset temperature
bldata.endset = Tendset;     % estimate of endset temperature
bldata.baselines = baselines;  % various other baselines

% select baseline via "type"
switch lower(bltype)
   case {'linrange', 'linear'}         % "linear" for backward compatibility
      blobj = bl_linrange;
   case {'linrange_sig', 'sigmoidal'}  % "sigmoidal" for backward compatibility
      blobj = bl_linrange_sig;
   case {'tote'} 
      blobj = bl_ToTe;
   case {'tote_sig'}
      blobj = bl_ToTe_sig;
   case {'simple'}
      blobj = bl_simple;
   case {'simple_sig'}
      blobj = bl_simple_sig;
   otherwise
      warning('Unknown baseline type requested: %s\nReturning linrange baseline', bltype)
      bltype = 'linrange';
      blobj = bl_linrange;
end

% finito
makeMessage('Reported default baseline type: "%s"\n', bltype);
return

end % of function