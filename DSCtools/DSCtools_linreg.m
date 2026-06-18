function reg = linreg(X, Y, reg)
% reg = linreg(X, Y, reg0)
% LINear REGression. 
%
% Assumes a linear model:  Yi = a + b*Xi + ei    with ei ~ N(0,s^2)
% and calculates the simple linear regression according to [1].
%
% INPUT:   X --> vector of x-values (scalars)
%          Y --> vector of y-values (scalars)
%       reg0 --> regression structure that shall be updated (may be empty)
%
% OUTPUT:  reg --> structure with following fields:
%           .a     --> estimate for a
%           .b     --> estimate for b
%           .s2    --> estimate for s^2
%           .Xmean --> mean value of Xi
%           .Ymean --> mean value of Yi
%           .Sxx   --> S_X^2
%           .Sxy   --> S_XY
%           .Syy   --> S_Y^2
%           .See   --> S_e^2
%
% [1] Jerome H. Klotz: "Updating Simple Linear Regression", 1995.
%     Statistica Sinica 5 (1995), 399-403
%
% Author: Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%

% error check
assert(length(X) == length(Y));

% if reg0 is empty, build reg from scratch
if (nargin < 3) || isempty(reg)
   X = reshape(X, [], 1);
   Y = reshape(Y, [], 1);
   n = length(X);
   Xmean = mean(X);
   Ymean = mean(Y);
   XminusXmean = X - Xmean;
   YminusYmean = Y - Ymean;
   Sxx = XminusXmean.' * XminusXmean;
   % reg.Syy = YminusYmean.' * YminusYmean;  % not needed
   Sxy = XminusXmean.' * YminusYmean;
   b   = Sxy / Sxx;
   a   = Ymean - b*Xmean;
   See = YminusYmean - b*XminusXmean;
   See = See.' * See;
   s2  = See / (n-2);
   % assign values
   reg.n     = n;
   reg.Xmean = Xmean;
   reg.Ymean = Ymean;
   reg.See   = See;
   reg.Sxx   = Sxx;
   reg.Sxy   = Sxy;
   reg.a     = a;
   reg.b     = b;
   reg.s2    = s2;
   return
end



% accessors
n      = reg.n;
bn     = reg.b;
%an    = reg.a;
Xmeann = reg.Xmean;
Ymeann = reg.Ymean;
Sxxn   = reg.Sxx;
Sxyn   = reg.Sxy;
Seen   = reg.See;



% if reg0 is given, update it
if isscalar(X)

   %% scalar update
   N = n + 1;
  
   XmeanN = Xmeann + (X - Xmeann)/N;
   YmeanN = Ymeann + (Y - Ymeann)/N;

   XminusXmean = X - Xmeann;
   YminusYmean = Y - Ymeann;

   SxxN = Sxxn + n/N * XminusXmean * XminusXmean;
   SxyN = Sxyn + n/N * XminusXmean * YminusYmean;
   SeeN = Seen + n/N * (YminusYmean - bn*XminusXmean)^2 * Sxxn / SxxN;
   
   bN = SxyN / SxxN;
   aN = YmeanN - bN*XmeanN;  % Formula (1)
   s2N = SeeN / (N-2);
   
else
   
   %% vector update
   X = reshape(X, [], 1);
   Y = reshape(Y, [], 1);

   nnew = length(X);
   N  = n+nnew;
   sqrtn = sqrt(n);
   sqrtN = sqrt(N);
   
   XmeanN = (n / N) * (Xmeann + nnew/n * mean(X));  % own derivation
   YmeanN = (n / N) * (Ymeann + nnew/n * mean(Y));  %
   
   XmeanSTAR = (Xmeann*sqrtn + XmeanN*sqrtN) / (sqrtn + sqrtN);  % below Formula (6)
   YmeanSTAR = (Ymeann*sqrtn + YmeanN*sqrtN) / (sqrtn + sqrtN);  % below Formula (7)
   
   XminusXmeanSTAR = (X - XmeanSTAR);               %
   SxxSTAR = XminusXmeanSTAR.' * XminusXmeanSTAR;   %
   SxxN    = Sxxn + SxxSTAR;                        % Formula (6)
   
   YminusYmeanSTAR = (Y - YmeanSTAR);               %
   SxySTAR = XminusXmeanSTAR.' * YminusYmeanSTAR;   %
   SxyN    = Sxyn + SxySTAR;                        % Formula (7)
   
   bN = SxyN / SxxN;                                % Formula (8)
   
   Sxn = sqrt(Sxxn);
   SxN = sqrt(SxxN);
   bSTAR = (bn * Sxn + bN*SxN) / (Sxn + SxN);          % Formula (10)
   
   SeeSTAR = YminusYmeanSTAR - bSTAR*XminusXmeanSTAR;  %
   SeeSTAR = SeeSTAR.' * SeeSTAR;                      % Formula (9)
   SeeN    = Seen + SeeSTAR;                           % Formula (8)
 
   aN  = YmeanN - bN*XmeanN;  % Formula (1)
   s2N = SeeN / (N-2);
   
end

%% store values
reg.n     = N;
reg.Xmean = XmeanN;
reg.Ymean = YmeanN;
reg.Sxx   = SxxN;
reg.Sxy   = SxyN;
reg.See   = SeeN;
reg.a     = aN;
reg.b     = bN;
reg.s2    = s2N;

end


