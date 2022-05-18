function sf = DSC204_getSensitivity(T)
   % sf = DSC204_getSensitivity(T)
   %
   % Retrieves the sensitivity factor for a certain temperature.
   %
   % The sensitivity factor sf is the reciprocal factor to transform
   % a milli-volt signal uV into a milli-watt power mW
   %             mW = uV / sf
   %

   % Measurements:
   % T_exp  = [ -48.3    -5.3    52.1    123.8  195.5    274.4    346.0    424.9    503.8    589.8   ];
   % K_calc = [   2.221   2.208   2.172    2.1    2.006    1.886    1.769    1.639    1.514    1.387 ];

   % set coefficients
   P0 = -5.3;
   P1 = 869.51691;
   P2 = 2.20833;
   P3 = -0.39066;
   P4 = -0.367;
   P5 = 1.3928;

   % calculate helper value z
   z = (T - P0) / P1;

   % calculate the factor
   sf = (P2 + P3*z + P4*z.^2 + P5*z.^3) .* exp(-z.^2);

end
