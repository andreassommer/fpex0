function cp = DSC204_cp_saphire_DIN11357(T, unit)
% cp = DSC204_cp_saphire_DIN11357(T)
%
% Delivers specific heat capacity of saphire for specified temperature according to DIN EN ISO 11357-4:2014-10
%
% INPUT:    T --> temperature
%        unit --> one of 'degC' or 'K'    (default: degC)
%
% OUTPUT:  cp --> literature value of saphire
%
% NOTE:   * This approximation is only valid in the interval 100K < T < 1200K
%
% Author: Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%


% unit not specified? then use degC
if (nargin < 2), unit = 'degC'; end


% coefficients in J/(g*K)
A = [  0.34407  ...
      -0.37623  ...
      -0.47824  ...
       0.54579  ...
       0.15393  ...
      -0.10023  ...
      -0.23778  ...
       0.26410  ...
      -0.21704  ...
       0.23260  ...
       1.12705 ];

    
% linear transformation of temperature
switch lower(unit)
   case {'degc','c','celsius'}; x = (T - 376.85) / 550;   % formula for degC
   case {'k','kelvin'};         x = (T - 650) / 550;      % formula for K
   otherwise
      error('Unknown unit: %s', unit)
end


% evaluate polynomial
cp = polyval(A, x);

end