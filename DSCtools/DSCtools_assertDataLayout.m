function layoutOK = DSCtools_assertDataLayout(DSCTOOLS_GLOBAL_SETTINGS)
% val = DSCtools_assertDataLayout(DSCTOOLS_GLOBAL_SETTINGS)
%
% Ensures that the specified data is in expected layout.
%
% INPUT:    DSCTOOLS_GLOBAL_SETTINGS --> DSCtools data structure as returned by DSCtools_readFile(s)
%
% OUTPUT:          val --> boolean indicating the data layout is as expected:
%                          {'Temp./DegC'  'Time/min'  'DSC/(uV/mg)'  'Sensit./(uV/mW)'}
%                          (where Deg denotes the degree symbol)
%

% compatibility with array call
if length(DSCTOOLS_GLOBAL_SETTINGS) > 1
   vals = arrayfun(@DSCtools_assertDataLayout, DSCTOOLS_GLOBAL_SETTINGS);
   layoutOK = all(vals);
   return
end

% assume the best
layoutOK = true;

% get encoding
enc = DSCTOOLS_GLOBAL_SETTINGS.fileEnc;

% prototype
expectedColumns = { sprintf('Temp./%cC', enc.degC),  'Time/min',  'DSC/(uV/mg)',  'Sensit./(uV/mW)' };

% check if colums are as expected
if ( length(expectedColumns) ~= length(DSCTOOLS_GLOBAL_SETTINGS.columnHeaders) )
   warning('FILE: %s  --- Column headers differ in length!', DSCTOOLS_GLOBAL_SETTINGS.fileSpec);
end

% check individual columns (let it explode, if too few)
for k = 1:length(expectedColumns)
   if not(strcmpi(expectedColumns{k}, DSCTOOLS_GLOBAL_SETTINGS.columnHeaders{k}))
      warning('FILE: %s  --- Column #%d, expected: %s, got: %s', DSCTOOLS_GLOBAL_SETTINGS.fileSpec, k, expectedColumns{k}, DSCTOOLS_GLOBAL_SETTINGS.columnHeaders{k});
      layoutOK = false;
   end
end



end