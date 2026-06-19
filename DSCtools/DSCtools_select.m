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
   %                         'id'       sub string (regexp) in ID field
   %                         'filespec' sub string (regexp) in file spec field
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
   [ ratemin, args] = olGetOption(args, 'ratemin', -inf());
   [ ratemax, args] = olGetOption(args, 'ratemax', +inf());
   [ massmin, args] = olGetOption(args, 'massmin', -inf());
   [ massmax, args] = olGetOption(args, 'massmax', +inf());
   [ regexp_id      , args] = olGetOption(args, 'id'       , '');
   [ regexp_filespec, args] = olGetOption(args, 'filespec' , '');

   % select by bounds
   rates  = [dscdata.rate];    dscdata = dscdata( ( rates >= ratemin) & ( rates <= ratemax) );
   masses = [dscdata.mass];    dscdata = dscdata( (masses >= massmin) & (masses <= massmax) );

   % regexp on ID
   if ~isempty(regexp_id)
      allIDs = {dscdata.ID};
      found = find_regexp(allIDs, regexp_id);
      dscdata = dscdata(found);
   end

   % regexp on filespec
   if ~isempty(regexp_filespec)
      allfiles = arrayfun(@(x) x.rawData.fileSpec, dscdata, 'uniformoutput', false);
      found = find_regexp(allfiles, regexp_filespec);
      dscdata = dscdata(found);
   end

   % finito -- unprocessed args left?
   olWarnIfNotEmpty(args)
   return

end


%% HELPER
function found = find_regexp(cellstring, expression)
   % delivers the indices of cell strings that contain the regular expression
   found = regexp(cellstring, expression);
   found = ~(cellfun(@isempty, found));   % = cellfun(@(x) ~isempty(x), found)
end