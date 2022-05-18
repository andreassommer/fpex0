function mW = DSC204_uV_to_mW(uV, T)
   % mW = DSC204_uV_to_mW(uV, T)
   %
   % Transforms a micro-volt signal at a certain temperature into a power in milli-watts.
   %
   % INPUT:    uV --> micro-volt signal
   %            T --> temperature at where uV was measured
   %
   % OUTPUT:   mW --> milli-watt power corresponding to the milli-watt signal
   %
   % This is based upon the sensitivity calibration as performed
   % by the software Proteus of the Netzsch DSC204 apparatus.
   %

   % retrieve the sensitivity factors
   sf = DSC204_getSensitivity(T);
   
   % apply them
   mW = uV ./ sf;

   % that's all
   
end


