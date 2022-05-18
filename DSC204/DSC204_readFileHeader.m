function [desc, enc, fileHeaderStrings] = DSC204_readFileHeader(fid, enc)
   % [desc, enc, fileHeaderStrings] = DSC204_readFileHeader(fid, enc)
   %
   % Reads the header (lines starting with # until the first empty line is found) 
   % of a DSC204 file opened in file handle fid.
   %
   % INPUT:   fid --> file id to DSC204 file
   %          enc --> encoding structure (defaults)
   %
   % OUTPUT:        desc --> descriptive structure (possibly array)
   %                 enc --> possibly updated encoding structre
   %   fileHeaderStrings --> cell array of strings containing the header lines
   %
   % Author: Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %
   
   % save current file position, rewind file
   filePosition = ftell(fid);
   frewind(fid);
   
   % use default encoding, if not given
   if (nargin<2) || isempty(enc)
      enc = DSC204_getEncodingDefaults();
   end
   
   % retrieve length of line ending markers (newline)
   [~, cutlen] = DSC204_getNewlineDesc(fid);
  
   % prepare output variable
   fileHeaderStrings = {};
   
   % scan the file until the first empty line or the ## tag is found
   while ~feof(fid)
      
      % read line, cut newline chars
      rawline = fgets(fid);
      rawlen  = length(rawline);
      line    = rawline(1:(rawlen-cutlen));
      linelen = rawlen - cutlen;
      
      % analyse line
      if (linelen==0)
         % empty line, continue reading
      elseif (linelen>=2) && strcmp(line(1:2),'##')
         % line begins with ## --> column headers, then rewind and stop reading!
         fseek(fid, -rawlen, 0); 
         break
      elseif strcmp(line(1),'#')
         % line starts with a # (marks a comment in header)
         fileHeaderStrings{end+1} = line;                                                    %#ok<AGROW>
      else
         warning('Data content found before header was read!')
         break
      end
      
   end
   
   
   % replace strange characters in file header
   fileHeaderStrings = DSC204_replaceChars(fileHeaderStrings, enc);
   
   
   % check if eof
   if feof(fid)
      warning('Unexpected end of file!')
   end
   
   
   % try to determine the actual delimiters
   % ======================================
   decimalDelimiter = '';
   fieldDelimiter   = '';
   for k = 1:length(fileHeaderStrings)
      % decimal delimiter
      if strncmp(line, '#DECIMAL', 8)
         decimalDelimiter = extractDelimiter(line);
      end
      if strncmp(line, '#SEPARATOR', 10)
         fieldDelimiter   = extractDelimiter(line);
      end
   end
   
   
   % load encoding defaults and update it
   % ====================================
   if ~isempty(decimalDelimiter)
      if ~strcmp(enc.decimalDelimiter, decimalDelimiter)
         fprintf('Updating decimal delimiter from "%s" to "%s"\n', enc.decimalDelimiter, decimalDelimiter);
         enc.decimalDelimiter = decimalDelimiter;
      end
   end
   if ~isempty(fieldDelimiter)
      if ~strcmp(enc.fieldDelimiter, fieldDelimiter)
         fprintf('Updating field delimiter from "%s" to "%s"\n', enc.fieldDelimiter, fieldDelimiter);
         enc.fieldDelimiter = fieldDelimiter;
      end
   end
   
   
   
   % detect the number of columns in file header
   % ===========================================
   headerColumnCounts = cellfun(@length,strfind(fileHeaderStrings,';')) + 1;
   headerColumnCount = getMostFrequentElement(headerColumnCounts);
   if ~(all(headerColumnCounts == headerColumnCount))
      disp('Detected header columns:')
      fprintf('%d ',headerColumnCounts);
      fprintf('\n')
      warning('Inconsistent number of fields in file header. Using: %d', headerColumnCount)
   end
   %
   % NOTE:  if we have more than 2 columns in the file header, 
   %        then we must store a structure array.
   
   
   % process the file header
   % =======================
   fieldNames    = cell(length(fileHeaderStrings),1);
   fieldContents = cell(length(fileHeaderStrings),headerColumnCount-1);
   for k = 1:length(fileHeaderStrings)
     
     % extract current line and assert that it is a string
     line = fileHeaderStrings{k};
     
     % dissect into strings (stored in cell array)
     content = textscan(line, '%s', 'delimiter', enc.fieldDelimiter,'WhiteSpace','');
     content = content{1};  % unwrap first cell layer
     
     % make field name: ensure it consists only of letters; and trim value string
     idx = isstrprop(content{1},'alphanum');
     fieldNames{k} = genvarname(content{1}(idx)); % ensure it's an valid variable name
     fieldContents(k,:) = content(2:end);
  
   end
   
   % store in desc structure
   fieldContents = cellfun(@strtrim, fieldContents, 'UniformOutput', false);
   desc = cell2struct(fieldContents, fieldNames);
   
 
   
   % finito: restore file position
   fseek(fid, filePosition, 'bof');
   return

   
   
   
   % ========================================
   
   % helper to translate delimiter
   % (not fully tested, as I don't know all possible strings)      
   function delim = extractDelimiter(str)
      if strfind(str, 'SEMICOLON'), delim = ';'; 
      elseif strfind(str, 'COMMA'), delim = ','; 
      elseif strfind(str, 'COLON'), delim = ':'; 
      elseif strfind(str, 'POINT'), delim = '.';
      elseif strfind(str, 'DOT')  , delim = '.'; 
      else
         delim = '';  % kind of error value
      end
   end
   
   % get maximum value of most frequent elements
   function val = getMostFrequentElement(array)
      array = sort(array, 'ascend');
      diffidx = find(diff([reshape(array, 1, []) inf]));
      ddiffidx = diff(diffidx);
      if isempty(ddiffidx)
         val = array(1);  % all values are the same
      else
         mostfrequents = array(diffidx(find(ddiffidx==max(ddiffidx))+1));
         val = max(mostfrequents);
      end
   end

   
end
      
