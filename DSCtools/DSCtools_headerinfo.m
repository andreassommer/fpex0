function [result, fnames] = DSCtools_headerinfo(DSCdata, fname, enc)
% [result, fnames] = DSCtools_headerinfo(DSCdata, fname, enc)
% 
% Retrieves and displays information from header
%
%  INPUT:   DSCdata --> DSC data structure as returned by DSCtools_readFile()
%             fname --> Header field to extract; if fname is [], a list of header fields is returned.
%               enc --> file encoding   [optional]
%
%  OUTPUT:   result --> field content or cell of char array with list of fields
%            fnames --> field names as valid matlab name 
%          
%
% Author:  Andreas Sommer, Jun2016
% code@andreas-sommer.eu
%

% quick check: must be scalar
if ~isscalar(DSCdata)
   error('This function only works on scalar DSCdata')
end

% initialize required variables
result = '';
fnames = '';

% use default encoding, if not given
if (nargin<3) || isempty(enc)
   enc = DSCtools_getEncodingDefaults();
end

% no fname given? 
if (nargin<2)
   fname = '';
end

% extract header
header = DSCdata(1).rawData.fileHeader;
n = length(header);

% get field names as strings from header
fstrings = cell(n, 1);
for i = 1:n  % walk through header lines
   tmp = extractBetween(header{i}, enc.fileHeaderSymbol, enc.headerDelimiter); % extract between # and :
   fstrings{i} = tmp{1};   % extractBetween returns a cell array of character vectors -- taking the first here
end

% special call: empty fname --> return list of field names
if isempty(fname)
   result = fstrings;
   fnames = matlab.lang.makeValidName(fstrings);
   return
end

% normal call: search for specified field name
foundIdx = 0;
for i = 1:n
   if strcmp(fname, fstrings{i})   % search for string match
      foundIdx = i;
      break
   end
end

% if exact string was not found, search for field name
if ~foundIdx
   fnames = matlab.lang.makeValidName(fstrings);
   for i = 1:n
      if strcmp(fname, fnames{i})  % search for field name match
         foundIdx = i;
         break
      end
   end
end

% return empty if not found
if ~foundIdx
   return
end

% get data
result = extractAfter(header{foundIdx}, enc.fieldDelimiter);  % get everything after the ; (field delimiter)
result = strtrim(result);   % remove whitespace

end % of function
