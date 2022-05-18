function sbc = DSC204_SeebeckCoefficients(T)
%
%
%

% Seebeck coefficients in uV/K
%
% Table 13 from DIN EN ISO 60584-1:2017-04, Chapter 6, page 15
% for thermocouple of type E
% 
sbc_data = [ ...
   -189.3442   27.7  ; ...  % Argon       TP (tripel point)
    -38.8344   54.0  ; ...  % Quecksilber TP
     +0.01     58.7  ; ...  % Wasser      TP
    +29.7646   61.4  ; ...  % Gallium     MP (melting point)
   +156.5985   71.6  ; ...  % Indium      FP (freezing point)
   +231.928    75.5  ; ...  % Zinn        FP
   +419.527    80.3  ; ...  % Zink        FP
   +660.323    80.1  ; ...  % Aluminum    FP
   +961.78     75.6  ];     % Silber      FP

% interpolate
sbc = interp1(sbc_data(:,1), sbc_data(:,2), T, 'spline');

% finito
end