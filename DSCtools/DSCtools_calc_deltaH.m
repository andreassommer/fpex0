function varargout = DSCtools_calc_deltaH(DSCsample, varargin)
% [deltaH, deltaHinfo] = DSCtools_calc_deltaH(DSCsample, varargin)
%
% Computes the latent heat of fusion of the DSCsample in multiple ways.
%
% INPUT:    DSCsample --> dsc data structure with added cp information
% 
% OUTPUT:      deltaH --> latent heat of fusion as integral over latent cp using default baseline
%          deltaHinfo --> calculate deltaH with other baselines for analysis
%
% NOTE:  If no output is requested, a message with deltaHinfo is displayed.
%
% Author: Andreas Sommer, Jul2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu

% arg checks
if (nargin < 1), error('Need input.'); end

% check if we have a cp field
if ~isfield(DSCsample, 'cp')
   error('DSC data does not have a cp field. Use DSCtools_addCP() first.')
end

% process args
args = varargin;
[doDisplay, args] = olGetOption(args, 'doDisplay', false);
olWarnIfNotEmpty(args);

% force display if no output variable is requested
doDisplay = doDisplay || (nargout == 0);

% output requested?
if (nargout > 0)
   deltaH_out     = zeros(size(DSCsample));
   deltaHinfo_out = cell(size(DSCsample));
end

% run through all dsc samples
for k = 1:numel(DSCsample)
   dk = DSCsample(k);

   % accessors to data
   T        = dk.cp.T;
   cpfun    = dk.cp.fun;
   cplatent = dk.cp.latentfun(T);

   % heat of fusion (integral over latentcp) -- for default baseline (used in in cplatent)
   deltaH = trapz(T, cplatent);

   % get list of available baselines
   baselines = dk.cp.bldata.baselines;
   blnames = fieldnames(baselines);

   % display all available baseline names
   if (doDisplay) && (k==1)
      makeMessage('List of baselines:\n')
      for i = 1:length(blnames)
         makeMessage('(%d) - %s\n', i, blnames{i});
      end
   end

   % walk through the baselines and compute the respective latent heat of fusion
   deltaHinfo = struct();
   for i = 1:length(blnames)
      baseline = baselines.(blnames{i});
      Tonset   = baseline.onset;
      Tendset  = baseline.endset;
      T_ToTe   = T(T>=Tonset & T<=Tendset);
      dHinfo.baseline = baseline;
      dHinfo.deltaH   = trapz(T_ToTe, cpfun(T_ToTe)-baseline.eval(T_ToTe));
      deltaHinfo.(blnames{i}) = dHinfo;  % store all in the deltaHinfo object
      if (doDisplay)
         if (i == 1)
            makeMessage('#%03d deltaH: ', k);
         end
         makeMessage('#raw', '  {%d} %6.2f   ', i, dHinfo.deltaH);
         if (i == length(blnames))
            makeMessage('#raw','  ID: %20s\n', dk.ID);
         end
      end
   end

   % store output (if requested)
   if (nargout > 0)
      deltaH_out(k) = deltaH;
      deltaHinfo_out{k} = deltaHinfo;
   end
end


% assemble output
if (nargout > 0)
   varargout{1} = deltaH_out;
   varargout{2} = [deltaHinfo_out{:}];  % change into struct array
end


end % of function
