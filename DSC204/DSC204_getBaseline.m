function [blfun, bldata] = DSC204_getBaseline(X, Y, type, blds)
% [blfun, bldata] = DSC204_getBaseline(X, Y, type, blds)
%
% Retrieves the baselevel function for specified DSC204data.
%
%  INPUT:   X --> vector of x-values (e.g. temperatures)
%           Y --> vector of y-values (e.g. cp-values)
%        type --> 'linear' or 'sigmoidal'                              (default: 'linear')
%        blds --> structure containing the BaseLine Detection Setting  (none for defaults)
%
%  OUTPUT:  blfun --> function handle of baselevel function
%          bldata --> structure containing information about the baseline
%
% Author:  Andreas Sommer, Apr2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%

% halt on warning?
global DSC204_getBaseline_haltOnWarning
if isempty(DSC204_getBaseline_haltOnWarning)
   DSC204_getBaseline_haltOnWarning = true;
end

% if no blds are specified, use defaults
if (nargin < 4),  blds = DSC204_getBaselineDetectionSettings(); end

% if type is not specified, set it to linear
if (nargin < 3),  type = 'linear'; end
 

% detect linear ranges
fprintf('Detecting linear ranges...')
[peakVal, peakPos] = max(Y);
initlen_L = floor(blds.L.initfraction * peakPos);
initlen_R = floor(blds.R.initfraction * (length(X)-peakPos));
detectionArgs_L = {blds.L.reldevA, blds.L.reldevB, blds.L.reldevS2, blds.L.absdevA, blds.L.absdevB, blds.L.absdevS2};
detectionArgs_R = {blds.R.reldevA, blds.R.reldevB, blds.R.reldevS2, blds.R.absdevA, blds.R.absdevB, blds.R.absdevS2};
[idx_L, reg_L] = DSC204_detectLinearRange(X,Y,'left' ,initlen_L,detectionArgs_L{:});
[idx_R, reg_R] = DSC204_detectLinearRange(X,Y,'right',initlen_R,detectionArgs_R{:});
fprintf('done.')

% get baselevel function
blfun = DSC204_getBaselinePrimitive(X, Y, idx_L, reg_L.a, reg_L.b, idx_R, reg_R.a, reg_R.b, type);


% small plausibility check:  slope of linear parts should be small;
maxslope = 0.1;
if (abs(reg_L.b) > maxslope) || (abs(reg_R.b) > maxslope)
   fprintf('\n');
   warning('Slope of linear part is large: left: %g, right: %g. There''s probably something wrong. Using proposed baseline instead!', reg_L.b, reg_R.b)
   oldblfun = blfun;
   if (abs(reg_L.b) > maxslope); reg_L.b = 0; reg_L.a = Y(1)  ; idx_L = 2          ; end
   if (abs(reg_R.b) > maxslope); reg_R.b = 0; reg_R.a = Y(end); idx_R = length(Y)-1; end
   blfun = DSC204_getBaselinePrimitive(X, Y, idx_L, reg_L.a, reg_L.b, idx_R, reg_R.a, reg_R.b, type);
   if DSC204_getBaseline_haltOnWarning
      figure(271); clf; plot(X,Y,X,oldblfun(X),X,blfun(X)); legend('data', 'detected baseline','proposed baseline', 'Location','best');
      disp('Entering debugger... Press F5 to proceed.');
      keyboard
   end
end

% DEBUG MODE:
%figure(271); clf; plot(X,Y,X,blfun(X)); legend('data', 'using baseline', 'Location','best'); disp('Entering debugger... Press F5 to proceed.'); keyboard


% if queried, build data structure
if (nargout >= 2)
   bldata.reg_L  = reg_L;     % regression information left
   bldata.reg_R  = reg_R;     % regression information right
   
   % classic onset/endset estimation: where the linear parts are left (poor estimation)
   % bldata.onset  = X(idx_L);
   % bldata.endset = X(idx_R);
   
   % better onset/endset estimation: point X where the data (X,Y) first falls below baseline (seen from peak maximum)
   % with fallback using the ind_L or idx_R
   bloffset = 0.02; % set to 0 for finding first below baseline
   idxOnset  = find( Y(  1:peakPos) - blfun(X(  1:peakPos)) < bloffset, 1, 'last' );
   idxEndset = find( Y(peakPos:end) - blfun(X(peakPos:end)) < bloffset, 1, 'first') + peakPos - 1;
   if isempty(idxOnset),  idxOnset  = idx_L; end
   if isempty(idxEndset), idxEndset = idx_R; end
   bldata.onset  = X(idxOnset );
   bldata.endset = X(idxEndset);
end

% finito
return




   


end