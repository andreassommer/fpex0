function layoutOK = DSC204_assertDataLayout(DSC204data)
% val = DSC204_assertDataLayout(DSC204data)
%
% Ensures that the specified data is in expected layout.
%
% INPUT:    DSC204data --> DSC204 data structure as returned by DSC204_readFile(s)
%
% OUTPUT:          val --> boolean indicating the data layout is as expected:
%                          {'Temp./DegC'  'Time/min'  'DSC/(uV/mg)'  'Sensit./(uV/mW)'}
%                          (where Deg denotes the degree symbol)
%

% compatibility with array call
if length(DSC204data) > 1
   vals = arrayfun(@DSC204_assertDataLayout, DSC204data);
   layoutOK = all(vals);
   return
end

% assume the best
layoutOK = true;

% get encoding
enc = DSC204data.fileEnc;

% prototype
expectedColumns = { sprintf('Temp./%cC', enc.degC),  'Time/min',  'DSC/(uV/mg)',  'Sensit./(uV/mW)' };

% check if colums are as expected
if ( length(expectedColumns) ~= length(DSC204data.columnHeaders) )
   warning('FILE: %s  --- Column headers differ in length!', DSC204data.fileSpec);
end

% check individual columns (let it explode, if too few)
for k = 1:length(expectedColumns)
   if not(strcmpi(expectedColumns{k}, DSC204data.columnHeaders{k}))
      warning('FILE: %s  --- Column #%d, expected: %s, got: %s', DSC204data.fileSpec, k, expectedColumns{k}, DSC204data.columnHeaders{k});
      layoutOK = false;
   end
end



end