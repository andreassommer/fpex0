function DSC204data = DSC204test_processing(DSC204data)
% DSC204data = DSC204_testprocessing(DSC204data)
%
% INPUT:   DSC204data --> data struct as returned by DSC204_readFile
%
% OUTPUT:  DSC204data --> with additional structure fields
%
% WARNING: this is just for testing!
%

% if DSC204data is an array, apply this function recursively 
if length(DSC204data) > 1
   DSC204data = arrayfun(@DSC204_testprocessing, DSC204data);
   return
end

% prepare persistent information
persistent DSCsaphire DSCsaphire_source

% defaults for saphire
saphire_filepath = DSC204_getDSCDataPath();
saphire_filename = 'ExpDat_16-398s-7.csv';
saphire_filespec = fullfile(saphire_filepath, saphire_filename);

% try to load the saphire data (one fixed rate at 10K/min)
if isempty(DSCsaphire) || ~strcmp(DSCsaphire_source, saphire_filespec)
   fprintf('Reading saphire measurements file %s ...', saphire_filespec);
   DSCsaphire = DSC204_readFile(saphire_filespec);
   DSCsaphire_source = saphire_filespec;
   fprintf('Done!\n')
end

% restrict temperature range (to ensure T is monotonically increasing)
Tmin = 50;
Tmax = 155;

% quick accessors to data
[tS, TS, uVS, sfS, mWS] = DSC204_quickAccessors(DSC204data, Tmin, Tmax);
tS = tS - tS(1);  % start with time 0

% quick accessors to saphire (Reference)
[tR, TR, uVR, sfR, mWR] = DSC204_quickAccessors(DSCsaphire, Tmin, Tmax);
%tR = tR - tR(1);

% make interpolations for saphire (at the points TS, i.e. temperature of the sample)
uVRi = interp1(TR, uVR, TS, 'linear');
%mWRi = interp1(TR, mWR, TS, 'linear');
%sfRi = interp1(TR, sfR, TS, 'linear');

% calculate cp values of s
cpR = DSC204_cp_saphire_DIN11357(TS);
cpS = DSC204_calc_cp_DIN11357(TS, 1, uVS, cpR, 1, uVRi, 0);
% NOTES: * both the uV-Signals are uV/mg, so already normalized to mass 1!
%        * we assume that both are "corrected" signale, i.e. the U0 is already substracted

% Normalize cp-value
cpS = cpS ./ DSC204data.Tinfo.Tstep;   % divide by Tdot --> normalizes to beta = 1 K/min
cpS = cpS *  DSCsaphire.Tinfo.Tstep;   % multiply by Tdot of saphire (now we have: beta_sample = beta_saphire)

% Normalize voltages
uVS_normalized = uVS ./ DSC204data.Tinfo.Tstep;

% store quick stuff in data structure
DSC204data.t  =  tS;  % time
DSC204data.T  =  TS;  % temperature
DSC204data.uV = uVS;  % DSC signal in micro-Volts
DSC204data.mW = mWS;  % heat flux in milli-Watts
DSC204data.sf = sfS;  % sensitivity factors in uV/mW
DSC204data.cp = cpS;  % apparent specific heat
DSC204data.uVn = uVS_normalized; 

end