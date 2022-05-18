function uV = DSC204_EMFwrt0(T)
%  uV = DSC204_EMFwrt0(T)
%
%  Delivers the thermoelectric voltage (EMF) in uV (micro-Volts) at temperature T (in degree Celsius),
%  at a reference temperature of 0 degree Celsius, for a type E thermocouple.
%
%  See DIN EN 60584-1:2014-07
%
%  INPUT:    T --> temperature in degree Celsius
% 
%  OUTPUT:  uV --> induced thermoelectric voltage in micro-Volts
%
%  Author: Andreas Sommer, May17
%  andreas.sommer@iwr.uni-heidelberg.de
%  email@andreas-sommer.eu
%
%

% coefficients if T in [-270, 0)
c_low = [
   +0.0d0
   +5.8665508708d+01
   +4.5410977124d-02
   -7.7998048686d-04
   -2.5800160843d-05
   -5.9452583057d-07
   -9.3214058667d-09
   -1.0287605534d-10
   -8.0370123621d-13
   -4.3979497391d-15
   -1.6414776355d-17
   -3.9673619516d-20
   -5.5827328721d-23
   -3.4657842013d-26 ];

% coefficients if T in [0, 1000]
c_high = [ 
   +0.0d0
   +5.8665508710d+01   
   +4.5032275582d-02
   +2.8908407212d-05
   -3.3056896652d-07
   +6.5024403270d-10
   -1.9197495504d-13
   -1.2536600497d-15
   +2.1489217569d-18
   -1.4388041782d-21
   +3.5960899481d-25 ];

n = length(c_high);
  
% horner schema  
uV = zeros(size(T));
idxHI = (T>=0);
idxLO = (T< 0);
uV(idxHI) = c_high(end);
uV(idxLO) = c_low(end);
for k = n:-1:1
   uV(idxHI) = T(idxHI) .* uV(idxHI)  + c_high(k);
   uV(idxLO) = T(idxLO) .* uV(idxLO)  + c_low(k);
end


% finito
end
