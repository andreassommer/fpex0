function DSCsample = DSC204_addCP(DSCsample, DSCreference, Tmin, Tmax, signalsource)
% DSCsample = DSC204_addCP(DSCsample, DSCreference, Tmin, Tmax, signalsource)
%
% Adds cp values to DSC204 data structure.
% 
% Calculation of cp is done with DSC204_calc_cp_DIN11357.
%
% INPUT:     DSCsample --> DSC204 data structure of sample(s) as returned by DSC204_readFile(s) - (e.g. pcm)
%         DSCreference --> DSC204 data structure of reference as returned by DSC204_readFile(s) - (e.g. saphire)
%                 Tmin --> lower temperature bound                             [default:  55 degC]
%                 Tmax --> upper temperature bound                             [default: 160 degC]
%         signalsource --> choose 'uV' for DSC voltage or 'mW' for DSC heat flux   [default: 'uV']
%
% OUTPUT:    DSC204data --> DSCsample structure with additional field cp.
%
% Author:  Andreas Sommer, Apr2017, Oct2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%
% TODO:   single DSCreference --> use automatic rescaling
%          multi DSCreference --> sort by DSC rate, and use respective reference measurements
%

% setup
debugMode = true;

% defaults
if (nargin < 5 || isempty(signalsource)); signalsource = 'uV'; end
if (nargin < 4 || isempty(Tmax))        ; Tmax = 160;          end
if (nargin < 3 || isempty(Tmin))        ; Tmin = 55;           end

% make function applicable for struct arrays
if length(DSCsample) > 1
   DSCsample = arrayfun(@(x) DSC204_addCP(x, DSCreference, Tmin, Tmax, signalsource), DSCsample);
   return
end

% display setup
if debugMode, fprintf('DSC204_addCP:  Using signal source "%s", Tmin=%g, Tmax=%g\n', signalsource, Tmin, Tmax); end

% retrieve the heating rates
betaS = DSCsample.rate;
betaR = [DSCreference.rate]; % possibly a vector

% reference rates not unique?
if not( length(betaR) == length(unique(betaR)) )
   fprintf('DSC204_addCP: Reference heat rates not unique! rate vector = [%s]\n', num2str(betaR));
   error('Reference heat rates not unique!')
end

% if multiple references are given, choose the one closest to the sample heat rate
if length(DSCreference) > 1
   if debugMode, fprintf('DSC204_addCP: Multiple references given; looking for heat rate beta=%g. ', betaS); end
   [~, refidx] = min(abs(betaR - betaS));
   DSCreference = DSCreference(refidx);
   betaR = DSCreference.rate;
   if debugMode, fprintf('Chosen reference rate: %g  (at index = %d)\n', betaR, refidx); end
end


% FROM HERE:  DSCreference is a SINGLE (ONE-ELEMENT) structure

% quick accessors
[~, TS, uVS, ~, mWS] = DSC204_quickAccessors(DSCsample);
[~, TR, uVR, ~, mWR] = DSC204_quickAccessors(DSCreference);

% message
fprintf('DSC204_addCP: Processing %-40s: ', DSCsample.fileSpec);

% assert data layout
if ( not(DSC204_assertDataLayout(DSCsample)) || not(DSC204_assertDataLayout(DSCreference)) )
   warning('Cannot ensure correct data layout. Proceed with fingers crossed!')
end


% select the signal type: µV (DSC raw signal) or mW (heat flux)
switch signalsource
   case 'uV'
      xS = uVS;
      xR = uVR;
   case 'mW'
      xS = mWS;
      xR = mWR;
   otherwise
      error('unknown signal source: %s', signalsource);
end


% determine Tmin and Tmax
Tmin = max( [min(TS), min(TR), Tmin] );
Tmax = min( [max(TS), max(TR), Tmax] );


% restrict temperatures and align everything at the temperature information of the sample measurement
% note: this is the temperature of the empty reference crucible, not of the sample itself
% get the signals of the reference corresponding to the temperatures of the sample.
% we use linear interpolation here
idxR = (TR >= Tmin  &  TR <= Tmax); 
idxS = (TS >= Tmin  &  TS <= Tmax);
TR = TR(idxR);
xR = xR(idxR);
xS = xS(idxS);
TS = TS(idxS);
xR = interp1(TR, xR, TS, 'linear', 'extrap');


% % issue a warning if sample temperatures leave the reference temperature range
% if (min(TS)<min(TR)  ||  max(TS)>max(TR))
%    warning('Warning: sample temperature leaves reference temperature range!')
% end

% masses
mS   = DSCsample.mass;     % mass of sample
mR   = DSCreference.mass;  % mass of reference

% measurements are normalized to uV/mg, so we recover the original signal by multiplying with mass
xS = mS * xS;
xR = mR * xR;

% from carefully looking at the measurement data, we see that the voltage signal is proportional
% to the heating rate, with proportionality constant 1.
% so we normalize both the sample and the reference signal to a heating rate of 1.0 K/min.
% NOTE: this does also not interfere if betaR==betaS
dscS = xS / betaS;
dscR = xR / betaR;
dsc0 = 0;  % disabled, since we already work with corrected signals

% now retrieve the reference cp values of saphire (unit: degC)
cpR = DSC204_cp_saphire_DIN11357(TS, 'degC');

% calculate cp of sample
cpS = DSC204_calc_cp_DIN11357(TS, mS, dscS, cpR, mR, dscR, dsc0);


% store cp values and associated temperatures
cpS = reshape(cpS, [], 1);
TS  = reshape(TS, [], 1); 
cp.values = cpS;
cp.T      = TS;


% store the piecewise polynomial with pchip interpolation
cp.pp     = interp1(TS, cpS, 'pchip', 'pp');
cp.fun    = @(T) cpevaler(T);
   % Helper function to evaluate cp
   function val = cpevaler(T)
      val = NaN(size(T));
      nanIdx = (T<Tmin) | (T>Tmax);
      val(~nanIdx) = ppval(cp.pp, T(~nanIdx));
      if any(nanIdx)
         warning('Warning: leaving interpolation temperature range [%g, %g]', Tmin, Tmax)
      end
   end

% get the baseline
[blfun, bldata] = DSC204_getBaseline(TS, cpS, 'linear');
cp.blfun  = blfun;
cp.bldata = bldata;

% build the latent cp function               DSC204_substractBaseline(X , Yin, blfun,clearzero,nonnegative, onset       , endset)
[latentCPvals, latentCPfun, latentCPfunpp] = DSC204_substractBaseline(TS, cpS, blfun, false   , true      , bldata.onset, bldata.endset);
%[latentCPvals, latentCPfun, latentCPfunpp] = DSC204_substractBaseline(TS, cpS, blfun, true    , true      , bldata.onset, bldata.endset);
cp.latentdata  = latentCPvals;
cp.latentfun   = latentCPfun;
cp.latentfunpp = latentCPfunpp;

% store it in DSC data
DSCsample.cp = cp;

% DEBUG: Show detected baseline
% xShowBaseline(DSCsample);

% message
fprintf('Done.\n')

end