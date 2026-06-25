function varargout = DSCtools_identifyColumns(columnHeaders, varargin)
% [col1, col2, ..., unit1, unit2, ...] = DSCtools_identifyColumns(columnHeaders, search1, search2, ...);
%
% Searches the column header for specified terms.
%
% INPUT:  columnHeader --> cell array of string or cellstr, linke in the columnHeaders field of DSC measurements
%              search1 --> 1st search string 
%              search2 --> 2nd search string  
%                      ...
%
% OUTPUT:  col1 --> column of 1st search string (empty [] if not found)
%          col2 --> colums of 2nd search string (empty [] if not found)
%               ...
%         unit1 --> unit in column1  (empty '' if not found)
%         unit2 --> unit in column2  (empty '' if not found)
%               ...
%
% Author:  Andreas Sommer, Jun2026
% code@andreas-sommer.eu

% import settings
global DSCTOOLS_GLOBAL_SETTINGS      %#ok<GVMIS>
warnIfColumnNotFound  = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'dataMAT_warn_if_column_not_found'  , true);
warnIfColumnNotUnique = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'dataMAT_warn_if_column_not_unique' , true);


% make this usable for multiple arguments
nqueries = length(varargin);
if nqueries > 1
   cols  = cell(nqueries, 1);
   units = cell(nqueries, 1);
   for i = 1:nqueries
      [ cols{i}, units{i} ] = DSCtools_identifyColumns(columnHeaders, varargin{i});
   end
   if (nargout <= nqueries)
      varargout = col;
   else % also units requested
      varargout = [ cols , units ];
   end
   return
end

%%% from here, we have a single keyword to search for
key = varargin{1};

% check in which columne the keyword can be found
idx = arrayfun(@(columnStr) contains(columnStr, key), columnHeaders);
idx = find(idx);  % get the linear indices

% if none found, warn
if isempty(idx) && warnIfColumnNotFound
   warning('Cannot find requested entry %s', key);
end

% if multiple found, warn
if numel(idx) > 1 && warnIfColumnNotUnique
   warning('Found multiple occations of %s: %s', key, sprintf('%d ', idx));
end

% return the indices
varargout{1} = idx;
if (nargout > 1)
   varargout{2} = getUnit(columnHeaders{idx});
end
return   


end % of function





%% HELPERS

function unit = getUnit(columnHeader)
   % Extract the unit
   slashpos = findCentralSlash(columnHeader);
   unit = columnHeader(slashpos+1:end);
   % remove a wrapping parentheses
   while hasSuperfluousOuterParens(unit)
      unit = unit(2:end-1);
   end
end


function idx = findCentralSlash(s)
% Find the central slash. We might have expressions like
% 'Gas Flow(protective)/(ml/min)'
% 'DSC/(mW/mg)'
% 'Temp./°C'
% 'Sensit./(uV/mW)'
   level = 0;
   idx = [];
   for k = 1:numel(s)
      if s(k) == '('
         level = level + 1;
      elseif s(k) == ')'
         level = max(level - 1, 0);
      elseif s(k) == '/' && level == 0
         idx = k;
         return
      end
   end
end


function tf = hasSuperfluousOuterParens(s)
% check if the outher parenthesis is superfluous
  if isempty(s) || s(1) ~= '(' || s(end) ~= ')'
     tf = false;
     return
  end
  level = 0;
  for k = 1:numel(s)
     if s(k) == '('
        level = level + 1;
     elseif s(k) == ')'
        level = level - 1;
        if level == 0 && k < numel(s)
           tf = false;
           return
        end
     end
  end
  tf = (level == 0);
end