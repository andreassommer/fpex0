function Tinfo = DSC204_extractRANGE(rangeString, enc)
   % Tinfo = DSC204_extractRANGE(rangeString, enc)
   %
   % Extracts start and end temperatures as well as the stepping.
   % If Tstart > Tend, then Tstep will become negative.
   %
   % INPUT:   rangeString --> string from field "RANGE" of the DSC204 file header
   %                          Typically, it looks like:  160°C/10,0(K/min)/25°C
   %                  enc --> encoding settings
   %
   %
   % OUTPUT:  Tinfo --> structure with fields:
   %             .Tstart --> initial temperature
   %              .Tstep --> stepping (negative for cooling)
   %               .Tend --> end temperature
   %          .Tstepunit --> detected unit of stepping (string as retrieved from rangeString)
   %              .Tunit --> detected temperature unit (string as retrieved from rangeString)
   %
   % Author: Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %
   
   % Range-String:  160°C/10,0(K/min)/25°C
   
   % first: trim whitespace
   rangeString = strtrim(rangeString);
   
   % find and remove the unit string (K/min)
   openParenIdx  = strfind(rangeString, '(');
   closeParenIdx = strfind(rangeString, ')');
   Tstepunit     = rangeString(openParenIdx+1:closeParenIdx-1);
   rangeString   = rangeString([1:openParenIdx-1  closeParenIdx+1:length(rangeString)]);
   
   % find and remove the temperature unit (heuristic: start from end until the first number is seen)
   lastDigitIdx = find(isstrprop(rangeString,'digit'),1,'last');
   Tunit        = rangeString(lastDigitIdx+1:end);
   rangeString  = strrep(rangeString, Tunit, '');
   
   % if decimal delimiter is not a point, replace it
   rangeString  = DSC204_replaceDecimalDelimiter(rangeString, enc);

   % now the string looks like: 160/10.0/25 and can be dissected
   Tcell  = textscan(rangeString, '%f', 'delimiter', '/');
   Tdata  = Tcell{1};  % unwrap first cell layer
   Tstart = Tdata(1);
   Tend   = Tdata(3);
   if Tstart > Tend
      Tstep = - Tdata(2);
   else 
      Tstep = Tdata(2);
   end
   
   % set the result structure
   Tinfo.Tstart    = Tstart;
   Tinfo.Tend      = Tend;
   Tinfo.Tunit     = Tunit;
   Tinfo.Tstep     = Tstep;
   Tinfo.Tstepunit = Tstepunit;

   
end