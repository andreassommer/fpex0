function dataMat = DSC204_readData(fid, enc, columnCount, emptyVal)
   % columnHeaders = DSC204_readData(fid, enc, columnCount)
   %
   % Imports data from DSC204 export file.
   %
   % INPUT:     fid --> file id to DSC204 file
   %            enc --> encoding structure
   %    columnCount --> number of columns to import
   %       emptyval --> value used for empty cells (default: NaN)
   %
   % OUTPUT:  dataMat --> cell array of strings
   %
   % Author: Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %

   % save current file position, rewind file
   filePosition = ftell(fid);
   frewind(fid);
   
   % retrieve length of line ending markers (newline)
   [~, cutlen] = DSC204_getNewlineDesc(fid);
      
   % prepare output variables
   replaceDecimalDelimimiter = ~strcmp(enc.decimalDelimiter, '.');
   data = {};
  
   % walk through file
   try
      linecount = 0;
      while ~feof(fid)
         % read line (this is 20% faster than using fgetl)
         rawline = fgets(fid);
         line    = rawline(1:length(rawline)-cutlen);
         linelen = length(line);
         linecount = linecount + 1;
         
         % analyse line
         if (linelen>0) 
            % skip commentary
            if (line(1)=='#'); continue; end
            % process data line
            if replaceDecimalDelimimiter
               %line = strrep(line,enc.decimalDelimiter,'.');
               line = DSC204_replaceDecimalDelimiter(line, enc);
            end
            % NOTE:  textscan has a problem if there is no value after the last delimiter.
            %        It does NOT substitute this with the "emptyvalue", but ignores it completely and reduces
            %        the number of read elements. So we have to add one.
            % read and interprete the line
            %datacells = textscan(line,'%f','Delimiter',enc.fieldDelimiter,'WhiteSpace'  ,'','EmptyValue',emptyVal);  % For Matlab2013
            %datacells = textscan(line,'%f','Delimiter',enc.fieldDelimiter,'TreatAsEmpty','','EmptyValue',emptyVal);  % For Matlab2017, breaks R2022a
            datacells = textscan(line,'%f','Delimiter',enc.fieldDelimiter,'EmptyValue',emptyVal);  % shoud work for all
            datarow = reshape(datacells{1}, 1, []); % unwrap first cell layers
            if length(datarow) < columnCount
               datarow(end+1) = emptyVal;                        %#ok<AGROW>
            end
            % store it
            data{end+1} = datarow;                               %#ok<AGROW>
         else
            % empty line: just skip
         end
      end
      
   catch err
      fprintf('Error while reading data line # %d: \n', linecount);
      fprintf('Line content: %s\n', line);
      fprintf('Catched error:\n');
      disp(err)
      fprintf('Stack trace:\n');
      disp(err.stack)
      rethrow(err)
   end
   
   % assemble data matrix
   dataMat = cell2mat(reshape(data,[],1));
 
   % finito: restore file position
   fseek(fid, filePosition, 'bof');
   return
   
end