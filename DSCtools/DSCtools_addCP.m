function DSCsample = DSCtools_addCP(DSCsample, DSCreference, Tmin, Tmax, signalsource, blds_preset, baselinetype)
% DSCsample = DSCtools_addCP(DSCsample, DSCreference, Tmin, Tmax, signalsource, blds_preset, baselinetype)
%
% Adds cp values to DSCtools data structure.
% 
% Calculation of cp is done with DSCtools_calc_cp_DIN11357.
%
% INPUT:     DSCsample --> DSCtools data structure of sample(s) as returned by DSCtools_readFile(s) - (e.g. pcm)
%         DSCreference --> DSCtools data structure of reference as returned by DSCtools_readFile(s) - (e.g. saphire)
%                          If DSCreference is empty, an artificial reference signal of permanently 1.0 is used.
%                 Tmin --> lower temperature bound                         [default taken from DSCTOOLS_GLPBAL_SETTINGS]
%                 Tmax --> upper temperature bound                         [default taken from DSCTOOLS_GLPBAL_SETTINGS]
%         signalsource --> choose 'uV' for DSC voltage or 'mW' for DSC heat flux   [default: 'uV']
%          blds_preset --> baseline detection settings preset        [optional, otherwise from DSCTOOLS_GLOBAL_SETTINGS]
%         baselinetype --> type of baseline                          [optional, overrides blds_preset]
%
% OUTPUT:    DSCsample --> DSCsample structure with additional field cp.
%
% NOTE:   DSCsample and DSCreference must be 0-corrected already.
%
% Author:  Andreas Sommer, Apr2017, Oct2017, Aug2022, Jun2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%
% TODO:   single DSCreference --> use automatic rescaling
%          multi DSCreference --> sort by DSC rate, and use respective reference measurements
%


% import settings
global DSCTOOLS_GLOBAL_SETTINGS                  %#ok<GVMIS>

debugMode   = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'addCP_debugMode'  );  % more verbosity?
scaleByMass = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'addCP_scaleByMass');  % multiply signaly by mass (signals are often normed to mass 1)
scaleToRate = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'addCP_scaleToRate');  % scale signals to heatrate (higher rates induce higher signals, see below)
TrangeWarn  = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'addCP_TrangeWarn' ); % warn if Tmin or Tmax exceed sample's or reference's temperature range

% process input args
if (nargin < 7 || isempty(baselinetype)); baselinetype = ''; end
if (nargin < 6 || isempty(blds_preset)) ;  blds_preset = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'blds_preset');               end
if (nargin < 5 || isempty(signalsource)); signalsource = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'addCP_defaultSignalSource'); end
if (nargin < 4 || isempty(Tmax))        ;         Tmax = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'addCP_defaultTmax');         end % previously: 160
if (nargin < 3 || isempty(Tmin))        ;         Tmin = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'addCP_defaultTmin');         end % previously:  55

% warning if output is not used
if (nargout==0)
   warning('Call to %s is in vain, because output is not stored by caller %s!', mfilename(), getCaller(0));
end
   
% make function applicable for struct arrays
if length(DSCsample) > 1
   DSCsample = arrayfun(@(x) DSCtools_addCP(x, DSCreference, Tmin, Tmax, signalsource, blds_preset, baselinetype), DSCsample);
   if (nargout==0)    % warning if output is not used
      warning('Call to %s is in vain, because output is not stored by caller %s!', mfilename(), getCaller(1));
   end
   return
end


% ---- from here, DSCsample is a SINGLE (1-ELEMENT) structure


% is a reference given?
if isempty(DSCreference)
   if debugMode, makeMessage('No reference signal given. Using artificial value 1.0\n'); end
   DSCreference = DSCsample;
   DSCreference.data.uV = ones(size(DSCreference.data.uV));
   DSCreference.data.mW = ones(size(DSCreference.data.mW));
   DSCreference.data.sf = ones(size(DSCreference.data.sf));
   DSCreference.mass    = 1.0;
   DSCreference.rate    = 1.0;
end

% display setup
if debugMode, makeMessage('Using signal source "%s", Tmin=%g, Tmax=%g\n', signalsource, Tmin, Tmax); end

% retrieve the heating rates
betaS = DSCsample.rate;
betaR = [DSCreference.rate]; % possibly a vector

% reference rates not unique?
if not( length(betaR) == length(unique(betaR)) )
   makeMessage('Reference heat rates not unique! rate vector = [%s]\n', num2str(betaR));
   error('Reference heat rates not unique!')
end

% if multiple references are given, choose the one closest to the sample heat rate
if length(DSCreference) > 1
   if debugMode, makeMessage('Multiple references given; looking for heat rate beta=%g\n', betaS); end
   [~, refidx] = min(abs(betaR - betaS));
   DSCreference = DSCreference(refidx);
   betaR = DSCreference.rate;
   if debugMode, makeMessage('Chosen reference rate: %g  (at index = %d)\n', betaR, refidx); end
end


% ---- from here, DSCreference is a SINGLE (1-ELEMENT) structure


makeMessage('Processing %s\n', DSCsample.ID);

% Variable naming: AB
%    A:    T -> Temperature   m -> mass   dsc -> dsc data
%    B:    S -> Sample    R -> Reference  
 

% quick accessors
TR = DSCreference.data.T;
TS = DSCsample.data.T;

% select the signal type: uV (DSC raw signal) or mW (heat flux)
switch signalsource
   case 'uV'
      dscS = DSCsample.data.uV;
      dscR = DSCreference.data.uV;
   case 'mW'
      dscS = DSCsample.data.mW;
      dscR = DSCreference.data.mW;
   otherwise
      error('unknown signal source: %s', signalsource);
end

% calculate minima and maxima of sample and reference temperatures
minTS = min(TS);   maxTS = max(TS);
minTR = min(TR);   maxTR = max(TR);

% issue a warning if Tmax or Tmin exceeds reference or sample temperature range
if ( Tmin < max(minTS,minTR)  ||  Tmax > min(maxTS,maxTR) )
   if debugMode
      makeMessage(' [ Tmin | minTS | minTR = %g | %g | %g  Tmax | maxTS | maxTR = %g | %g | %g ]\n', ...
                      Tmin, minTS, minTR, Tmax, maxTS, maxTR);
   end
   if TrangeWarn
      warning('Specified Tmin or Tmax exceeds sample or reference temperature range. Will be adjusted!');
   end
end


% determine Tmin and Tmax
Tmin = max( [minTS, minTR, Tmin] );
Tmax = min( [maxTS, maxTR, Tmax] );

% restrict temperatures and align everything at the temperature information of the sample measurement
% note: this is the temperature of the empty reference crucible, not of the sample itself
% get the signals of the reference corresponding to the temperatures of the sample.
% we use linear interpolation and extrapolation here
idxR = (TR >= Tmin  &  TR <= Tmax); 
idxS = (TS >= Tmin  &  TS <= Tmax);
dscR = dscR(idxR);
dscS = dscS(idxS);
TR   = TR(idxR);
TS   = TS(idxS);
dscR = interp1(TR, dscR, TS, 'linear', 'extrap');

% masses
mS = DSCsample.mass;     % mass of sample
mR = DSCreference.mass;  % mass of reference

% measurements are normalized to uV/mg, so we recover the original signal by multiplying with mass
if scaleByMass
   dscS = mS * dscS;
   dscR = mR * dscR;
end
dsc0 = 0;   % disabled, since we already work with corrected signals

% from carefully looking at the measurement data, we see that the voltage signal is proportional
% to the heating rate, with proportionality constant approximately 1.
% so we normalize both the sample and the reference signal to a heating rate of 1.0 K/min.
% NOTE: this does also not interfere if betaR==betaS
if scaleToRate
   dscS = dscS / betaS;
   dscR = dscR / betaR;
end

% now retrieve the reference cp values of saphire (unit: degC)
cpR = DSCtools_cp_saphire_DIN11357(TS, 'degC');

% calculate cp of sample
cpS = DSCtools_calc_cp_DIN11357(TS, mS, dscS, cpR, mR, dscR, dsc0);

% store cp values and associated temperatures
cpS = reshape(cpS, [], 1);
TS  = reshape(TS, [], 1); 
cp.values = cpS;
cp.T      = TS;

% store the piecewise polynomial with pchip interpolation
cp.pp     = interp1(TS, cpS, 'pchip', 'pp');
cp.fun    = @(T) DSCtools_addCP_cpevaler(T, cp, Tmin, Tmax);

% get the baseline
[baseline, bldata] = DSCtools_getBaseline(TS, cpS, baselinetype, blds_preset);
cp.baseline = baseline;
cp.bldata   = bldata;

% build the latent cp function               DSCtools_subtractBaseline(X , Yin, blobj   , clearzero, nonnegative, onset, endset)
[latentCPvals, latentCPfun, latentCPfunpp] = DSCtools_subtractBaseline(TS, cpS, baseline, false    , true      );
cp.latentdata  = latentCPvals;
cp.latentfun   = latentCPfun;
cp.latentfunpp = latentCPfunpp;

% calculate latent heat of fusion (integral over latentcp, "Delta H")
T_ToTe   = cp.T(cp.T>=baseline.onset & cp.T<=baseline.endset);
cplatent = cp.latentfun(T_ToTe);
deltaH = trapz(T_ToTe, cplatent); 
cp.deltaH = deltaH;   % store it

% store in DSC data struct
DSCsample.cp = cp;

% DEBUG: Show detected baseline
% xShowBaseline(DSCsample);

% message
makeMessage('Done processing %s.\n', DSCsample.ID);

% finito
return

end % of function

