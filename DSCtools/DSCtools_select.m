function [dscdata, idx] = DSCtools_select(dscdata, varargin)
   % function [dscdata, idx] = DSCtools_select(dscdata, key-value-pair*)
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
   %                 idx --> indices of the selected subset in the input dscdata struct array
   %
   %
   % Author: Andreas Sommer  --  Mar2017, Aug2022, Jun2026
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %

   % quick return on empty data
   if isempty(dscdata), idx = []; return; end

   % process arguments
   args = varargin;
   [ ratemin        , args] = olGetOption(args, 'ratemin' , -inf());
   [ ratemax        , args] = olGetOption(args, 'ratemax' , +inf());
   [ massmin        , args] = olGetOption(args, 'massmin' , -inf());
   [ massmax        , args] = olGetOption(args, 'massmax' , +inf());
   [ regexp_id      , args] = olGetOption(args, 'id'      , '');
   [ regexp_filespec, args] = olGetOption(args, 'filespec', '');
   olWarnIfNotEmpty(args, true); % unprocessed args left?

   % select by rate
   rates    = [dscdata.rate];
   idx_rate = ( rates >= ratemin) & ( rates <= ratemax) ;
         
   % select by mass
   masses   = [dscdata.mass];
   idx_mass = (masses >= massmin) & (masses <= massmax) ;
   
   % regexp on ID
   idx_ID = false(size(dscdata));
   if ~isempty(regexp_id)
      allIDs = {dscdata.ID};
      linidx_ID = find_regexp(allIDs, regexp_id);
      idx_ID(linidx_ID) = true;
   else
      idx_ID = true;
   end

   % regexp on filespec
   idx_file = false(size(dscdata));
   if ~isempty(regexp_filespec)
      allfiles = arrayfun(@(x) x.rawData.fileSpec, dscdata, 'uniformoutput', false);
      linidx_file = find_regexp(allfiles, regexp_filespec);
      idx_file(linidx_file) = true;
   else
      idx_file = true;
   end

   % combine the respective indices and get the data
   idx = reshape(idx_rate, [], 1) ...
       & reshape(idx_mass, [], 1) ...
       & reshape(idx_ID  , [], 1) ...
       & reshape(idx_file, [], 1) ;
   dscdata  = dscdata( idx );

   % finito
   return

end


%% HELPER
function found = find_regexp(cellstring, expression)
   % delivers the indices of cell strings that contain the regular expression
   found = regexp(cellstring, expression);
   found = ~(cellfun(@isempty, found));   % = cellfun(@(x) ~isempty(x), found)
end