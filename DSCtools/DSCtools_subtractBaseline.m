function [Yvals, Yfun, Yfunpp] = DSCtools_subtractBaseline(X, Yin, baseline, clearzero, nonnegative, onset, endset)
   % [Yvals, Yfun, Yfunpp] = DSCtools_subtractBaseline(X, Yin, baseline, clearzero, nonnegative, onset, endset)
   % 
   % Subtracts the baseline from given points.
   %
   % INPUT:     X --> x values (e.g. vector temperatures)
   %            Y --> y values or function (e.g. vector of cp values, or function cp(T))
   %     baseline --> baseline object
   %    clearzero --> flag indicating to clear zeros (see DSCtools_clearZeroFromMax) [default: true]
   %  nonnegative --> flag indicating to ensure nonnegativity                        [default: true]
   %        onset --> onset value (zero values are put below/left of this x value)   [optional, overrides baseline.onset]
   %       endset --> endset value (zero values are put above/right of this x value) [optional, overrides baseline.endsset]
   %
   % OUTPUT   Yvals --> processed y values
   %           Yfun --> function of processed values
   %         YfunPP --> piecewise polynomial representation of Yfun
   %
   % Author:  Andreas Sommer, Apr2017, Jul2026
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   
   % defaults
   if (nargin < 4), clearzero   = true; end
   if (nargin < 5), nonnegative = true; end
   if (nargin < 6), onset  = baseline.onset ; end
   if (nargin < 7), endset = baseline.endset; end

   % if Y is a function, evaluate it at X, otherwise ensure correct dimensions
   if isnumeric(Yin)
      assert(all(size(X)==size(Yin)));
      Yvals = Yin;
   elseif isa(Yin,'function_handle')
      Yvals = Yin(X);
   else
      error('Cannot process Y of class %s.', class(Yin))
   end
   
   % substract baseline from Y data
   Yvals = Yvals - baseline.eval(X);
   
   % make zeros outside the interval [onset, endset]
   Yvals(X<onset) = 0;
   Yvals(X>endset) = 0;
   
   % nonnegativity
   if nonnegative
      Yvals(Yvals<0) = 0;
   end
   
   % clear zero
   if clearzero
      Yvals = DSCtools_clearZeroFromMax(Yvals);
   end
   
   % for the interpolation function, add some zero space left and right
   addlen = 5;
   XX = [X(1) + (-addlen:-1)  ,  reshape(    X, 1, [])  ,  X(end) + (1:addlen)];
   YY = [zeros(1, addlen)     ,  reshape(Yvals, 1, [])  ,  zeros(1, addlen)   ];
   
   % build function from values
   Yfunpp = interp1(XX,YY,'pchip','pp');
   Yfun = @(T) ppval(Yfunpp, T);
      
   % finito
   return
      
end