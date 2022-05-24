function blfun = DSC204_getBaselinePrimitive(X, Y, index_l, icept_l, slope_l, index_r, icept_r, slope_r, type, res)
% blfun = DSC204_getBaselinePrimitive(X, Y, index_l, icept_l, slope_l, index_r, icept_r, slope_r, type, res)
%
% Retrieves the baselevel function for specified data.
%
%  INPUT:    X  -->  x-values (temperatures, time, ...)
%            Y  -->  y-values (voltage, heat flux, ....)
%      index_l  -->  index where left linear part is left ("onset")
%      icept_l  -->  y-intercept of left linear part
%      slope_l  -->  slope of left linear part
%      index_r  -->  index where right linear part is left ("onset")
%      icept_r  -->  y-intercept of right linear part
%      slope_r  -->  slope of right linear part
%         type  -->  'linear' or 'sigmoidal'
%          res  -->  number of support points for sigmoidal (default: 100)
%
%  OUTPUT:  blfun  -->  function handle of baselevel function
%
% Author:  Andreas Sommer, Apr2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%

% if res is not specified, set default value
if (nargin < 10)
   res = 1000;
end
   

% calculate values at "onset" and "offset"
yval_l = icept_l + X(index_l) * slope_l;
yval_r = icept_r + X(index_r) * slope_r;
yval_delta = yval_r - yval_l;

% calculate linear base level function
lin_bl_xvals = [        X(1)        , X(index_l), X(index_r),          X(end)         ];
lin_bl_yvals = [icept_l+X(1)*slope_l,   yval_l  ,   yval_r  ,  icept_r+X(end)*slope_r ];
lin_bl_pp  = interp1(lin_bl_xvals, lin_bl_yvals, 'linear', 'pp');
lin_bl_fun = @(x) ppval(lin_bl_pp, x);

   
% select output as requested
switch lower(type)
   
   case 'linear'
      % just return the linear base level function
      blfun = lin_bl_fun;
   
   case {'sigmoid','sigmoidal'}
      % sigmoidal interpolation as in Fig. 3 of DIN 51007
      
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
      sig_x = [            X(1)              X(index_l-1)   sig_nodes                X(index_r+1)              X(end) ];
      sig_y = [ lin_bl_fun(X(1))  lin_bl_fun(X(index_l-1))  sig_nodevals  lin_bl_fun(X(index_r+1))  lin_bl_fun(X(end)) ]; 
      sig_pp = pchip(sig_x, sig_y);
      
      % base level function
      blfun = @(x) ppval(sig_pp, x);

   otherwise
      error('Unknown baseline type: %s', type)
end



% finito
return





end