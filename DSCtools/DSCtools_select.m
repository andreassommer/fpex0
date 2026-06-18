function dscdata = DSCtools_select(dscdata, varargin)
   % function dscdata = DSCtools_select(dscdata, key-value-pair*)
   %
   % Selects a subset of the dsc measurements.
   %
   % INPUT:      dscdata --> Struct array with contents as read by DSCtools_readFile()
   %      key-value-pair --> possible keys:
   %                         'ratemin'  lower bound for heat rate 
   %                         'ratemax'  upper bound for heat rate
   %                         'massmin'  lower bound for sample mass
   %                         'massmax'  upper bound for sample mass
   %
   % OUTPUT:     dscdata --> Struct array with contents as read by DSCtools_readFile.
   %                         See there for documentation.
   %
   %
   % Author: Andreas Sommer  --  Mar2017, Aug2022, Jun2026
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %

   % store the arguments
   args = varargin;

   % process the bound fields (if specified)
   [ratemin, args] = olGetOption(args, 'ratemin', -inf());
   [ratemax, args] = olGetOption(args, 'ratemax', +inf());
   [massmin, args] = olGetOption(args, 'massmin', -inf());
   [massmax, args] = olGetOption(args, 'massmax', +inf());

   % select by bounds
   rates  = [dscdata.rate];    dscdata = dscdata( ( rates >= ratemin) & ( rates <= ratemax) );
   masses = [dscdata.mass];    dscdata = dscdata( (masses >= massmin) & (masses <= massmax) );

   % finito -- unprocessed args left?
   olWarnIfNotEmpty(args)
   return

end