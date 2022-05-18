function Tout = DSC204_Tcorr(Tin)
   % Tout = DSC_Tcorr(Tin)
   %
   % Temperature correction as performed by Proteus software of the DSC204 apparatus.
   %
   % This is the correction due to the temperature calibration
   % ("Temperaturkalibrierung").

   % Values exported from the Proteus software
   % T_nom = [-64.5 156.6 231.9 271.4 419.5];    % nominal values from literature
   % T_exp = [-64.3 156.3 231.6 271.1 418.9];    % experimentally determined values
   % T_nom = [-64.5 156.6 231.8 271.4 419.5];    % unclear
   % T_cor = [-0.207348699 0.105459456 0.234920425 0.307711128 0.608812244];
   %          correction values for the above temperatures after fitting
       
   % Set coefficients
   B0 = -126.7;
   B1 =  132.1;
   B2 =  103.8;
   
   % Evaluate
   Tout = 1.0d-3 .* B0  +  1.0d-5 .* B1 .* Tin  +  1.0d-8 * B2 * Tin.^2;

end