function cp = DSC204_calc_cp_by_integration(DSC204meas)
   % Calculates apparent cp-values from 
   %
   % INPUT:    DSC204meas --> DSC204 data struct array as by DSC204_readfiles();
   %
   % OUTPUT:   cp --> apparent cp values
   %                  if input is a struct array, the result is a cell array.
   %
   % Seriously, we could as well just multiply by the heating/cooling rate!
   %
   
   warning('I still believe, this is a wrong approach!')
   
   % temperature and time range restriction
   % Tmin =  70; Tmax = 150;
   Tmin = -inf; Tmax = +inf;
   
   % quick accessors
   T  = DSC204meas.data(:,1);   % temperature assumed in degC
   t  = DSC204meas.data(:,2);   % time        assumed in min
   uV = DSC204meas.data(:,3);   % DSC-signal  assumed in micro-Volts
   sf = DSC204meas.data(:,4);   % Sensitivity assumed in milli-Watts / micro-Volt

   % restrict ranges
   idx = (T>=Tmin) & (T<=Tmax); T = T(idx); t = t(idx); uV = uV(idx); sf = sf(idx);
   
   % derive milli-Watt signals from micro-Vols DSC-signal using the sensitivity
   Qdot = uV ./ sf;             % Qdot = power(t) = DSCuV(t)/sensitivity(t)        [Qdot] = W = J/s

   % as Qdot is in J/s, we must transform time from minutes to seconds
   t = t * 60;
   
   % cumulative integration w.r.t. time
   Q = cumtrapz(t, Qdot);
   
   % take derivative w.r.t. temperature
   cp = gradient(Q, T);
   
end