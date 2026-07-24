function blobj = DSCtools_getBaselinePrimitive(X, Y, index_l, icept_l, slope_l, index_r, icept_r, slope_r, type, res)
% blobj = DSCtools_getBaselinePrimitive(X, Y, index_l, icept_l, slope_l, index_r, icept_r, slope_r, type, res)
%
% Retrieves the baselevel function for specified data.
%
%  INPUT:    X  -->  x-values (temperatures, time, ...)
%            Y  -->  y-values (voltage, heat flux, ....)
%      index_l  -->  index where left linear part is left ("onset")
%      icept_l  -->  y-intercept of left linear part
%      slope_l  -->  slope of left linear part
%      index_r  -->  index where right linear part starts ("endset")
%      icept_r  -->  y-intercept of right linear part
%      slope_r  -->  slope of right linear part
%         type  -->  'linear' or 'sigmoidal'
%          res  -->  number of support points for sigmoidal (default: taken from DSCTOOLS_GLOBAL_SETTINGS)
%
%  OUTPUT:  blobj  -->  baseline object of type DSCtools_Baseline
%
% Author:  Andreas Sommer, Apr2017, Jul2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%

% if res is not specified, set default value
global DSCTOOLS_GLOBAL_SETTINGS   %#ok<GVMIS>
if (nargin < 10)
   res = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'getBaselinePrimitive_resolution');
end
   

% calculate values at "onset" and "offset"
yval_l = icept_l + X(index_l) * slope_l;
yval_r = icept_r + X(index_r) * slope_r;
yval_delta = yval_r - yval_l;

% calculate linear base level function
lin_bl_xvals = [         X(1)          ;  X(index_l)  ;  X(index_r)  ;          X(end)         ];
lin_bl_yvals = [ icept_l+X(1)*slope_l  ;    yval_l    ;    yval_r    ;  icept_r+X(end)*slope_r ];
lin_bl_xvals = ensure_unique(lin_bl_xvals, 1);  % if index_l == 1 or index_r==end, we need to have unique values
lin_bl_pp  = pwlinpp([lin_bl_xvals , lin_bl_yvals]);
lin_bl_fun = @(x) ppval(lin_bl_pp, x);

% lin_bl_pp = interp1(lin_bl_xvals, lin_bl_yvals, 'linear', 'pp');  % old call, much slower

   
% select output as requested
switch lower(type)
   
   case 'linear'
      % just return the linear base level function
      blfun = lin_bl_fun;
   
   case {'sigmoid','sigmoidal'}
      % sigmoidal interpolation as in Fig. 3 of DIN 51007
      type = 'sigmoidal';  % unify for baseline object below

      % subtract baseline and integrate peak part "in between"
      Xmid = X(index_l:index_r);
      Ymid = Y(index_l:index_r) - lin_bl_fun(Xmid);
      Ymid(Ymid<0) = 0;
      cumarea = cumtrapz(Xmid,Ymid);  % cumulative integral (area)
      sigmoidal = cumarea / max(cumarea) * yval_delta + yval_l;  % should be at end, but who knows...
      
      % interpolate integral at support points (res = #intervals)
      step = (Xmid(end)-Xmid(1)) / res;
      sig_nodes = [Xmid(1)  (Xmid(1)+step):step:(Xmid(end)-step)  Xmid(end)];
      sig_nodevals = interp1(Xmid, sigmoidal, sig_nodes, 'linear');
      
      % generate baseline function (piecewise cubic hermite interpolant)
      if (index_l) <= 1        , index_l = 2          ; end
      if (index_r) >= length(X), index_r = length(X)-1; end
      sig_x = [            X(1)              X(index_l-1)   sig_nodes                X(index_r+1)              X(end) ];
      sig_y = [ lin_bl_fun(X(1))  lin_bl_fun(X(index_l-1))  sig_nodevals  lin_bl_fun(X(index_r+1))  lin_bl_fun(X(end))]; 
      sig_x = ensure_unique(sig_x); % pchip requires unique x points
      sig_pp = pchip(sig_x, sig_y);
      
      % base level function
      blfun = @(x) ppval(sig_pp, x);

   otherwise
      error('Unknown baseline type: %s', type)
end

% assemble the object
blobj = DSCtools_Baseline(blfun, X(index_l), X(index_r), type);

% finito
return





end