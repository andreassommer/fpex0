function [stopidx,  reg, initreg] = DSC204_detectLinearRange(X,Y,side,initlen,reldevAmax,reldevBmax,reldevS2max,absdevAmax,absdevBmax,absdevS2max)
% [stopidx, reg, initreg] = DSC204_detectLinearRange(X,Y,side,initlen,reldevA,reldevB,reldevS2,absdevA,absdevB,absdevS2)
% Detect the stop position of the linear range.
%
% INPUT:   X --> x-values
%          Y --> y-values
%       side --> 'L' or 'R' for "from left" or "from right"
%    initlen --> length (in samples) to calculate the initial standard deviation
%    reldevA --> acceptable relative deviation of initial y-axis intercept
%    reldevB --> acceptable relative deviation of initial slope
%   reldevS2 --> acceptable relative deviation of initial squared standard deviation
%    absdevA --> acceptable absolute deviation of initial y-axis intercept
%    absdevB --> acceptable absolute deviation of initial slope
%   absdevS2 --> acceptable absolute deviation of initial squared standard deviation
%
% OUTPUT: stopidx --> position in X, where the linear range is estimated to be left
%             reg --> final regression structure as returned by DSC204_linreg();
%         initreg --> initial regression structure as returned by DSC204_linreg();
%
% NOTES: * Set the reldev* to inf to disable
%        * Set the absdev* to  -1 to disable
%        * The reldev* is only be checked if value exceeds absdev*
%
% Author:  Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%


% reshapeing
X = reshape(X,[],1);
Y = reshape(Y,[],1);
assert(length(X)==length(Y));

% length 
lenX = length(X);

% testweise kann man auch mal die daten hier rumdrehen, falls side=='R'
% if strcmpi(side,'R')
%    side = 'L';
%    X=X(end:-1:1);
%    Y=Y(end:-1:1);
% end


% choose direction
switch upper(side(1))
   case 'L'
      fromLeft = 1;
   case 'R'
      fromLeft = 0;
   otherwise 
      error('Unknown side: %s. Use ''L'' or ''R''', side);
end


% select indices
if fromLeft
   idx = 1:initlen;
   runidx = initlen+1:lenX;
else
   idx = lenX:-1:lenX-(initlen-1);
   runidx = (lenX-(initlen-1)-1):-1:1;
end


% initial linear regression
initreg = DSC204_linreg(X(idx),Y(idx));
reg = initreg;

% walk through data: blockwise? maybe later.
for pos = runidx
   reg = DSC204_linreg(X(pos),Y(pos),reg);  % update regression
   absdevA  = abs(reg.a -initreg.a);   reldevA  = abs(absdevA  / initreg.a );
   absdevB  = abs(reg.b -initreg.b);   reldevB  = abs(absdevB  / initreg.b );
   absdevS2 = abs(reg.s2-initreg.s2);  reldevS2 = abs(absdevS2 / initreg.s2);
   % DEBUG
%    fprintf('POS = %5.d    a: %7.4f, %7.4f, %7.4f, %4.2f        b: %7.4f, %7.4f, %7.4f, %4.2f        s: %7.4g, %7.4g, %7.4f, %4.2f   \n', ...
%             pos,          reg.a, initreg.a, absdevA, reldevA,  reg.b, initreg.b, absdevB, reldevB,  reg.s2, initreg.s2, absdevS2, reldevS2);
   % /DEBUG
   
   % check deviation in A (y-axis intercept)
   if (absdevA  > absdevAmax ) && (reldevA  > reldevAmax ); break; end
   if (absdevB  > absdevBmax ) && (reldevB  > reldevBmax ); break; end
   if (absdevS2 > absdevS2max) && (reldevS2 > reldevS2max); break; end
   
end

% results
stopidx = pos;


% finito
end
