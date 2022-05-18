function DSC204test_analyseData_oKmK(fileSpec)

% config
figNum = 400;
doExport = false;


% fileSpec angegeben?
if (nargin<1)
   % Datenauswahl (nur, falls fileSpec leer ist!)
   filePath   = '~/Projekte/modELTES/DSC204_F1_Phoenix_Messungen/alles';
   fileSeries = '407';
   fileDir    = 'H';     % 'K' for cooling
   fileSpeed  = '1,25';
   fileProto  = sprintf('%s/ExpDat_16-%s-3_xxxKorr_%sKmin_%s.csv', filePath, fileSeries, fileSpeed, fileDir);
else
   % Ersetze 'ohne' und 'mit' durch 'xxx'
   fileProto  = fileSpec;
   fileProto  = strrep(fileProto, 'ohne', 'xxx');
   fileProto  = strrep(fileProto, 'mit', 'xxx');
end


% Lade Daten
fileSpec_oK = strrep(fileProto, 'xxx', 'ohne');
fileSpec_mK = strrep(fileProto, 'xxx', 'mit');
oK = DSC204_readFile(fileSpec_oK);
mK = DSC204_readFile(fileSpec_mK);

% Schnellzugriffe  oK = ohneKorrektur;  mK = mitKorrektur
oK.T  = oK.data(:,1);
oK.t  = oK.data(:,2);
oK.uV = oK.data(:,3);
oK.sf = oK.data(:,4);
mK.T  = mK.data(:,1);     % temperature in degC
mK.t  = mK.data(:,2);     % time in minutes
mK.uV = mK.data(:,3);     % DSC signal
mK.sf = mK.data(:,4);     % sensitivity factor
t = oK.t;
T = oK.T;

% check: t und T müssen übereinstimmen
assert(all(oK.T == mK.T));
assert(all(oK.t == mK.t));



%% KORREKTUREN
figNum = figNum + 1; figure(figNum); clf(figNum); hold(gca,'on');
plotyy(t, [oK.uV  mK.uV], t, oK.uV-mK.uV ); 
timelabel();
legend('oK', 'mK', 'diff', 'location', 'best')

%% SENSITIVITAETEN

% Vergleiche Sensitivitäten:  selbst berechnet und jeweils aus den Daten
oK.sf_calc = DSC204_getSensitivity(oK.T);
mK.sf_calc = DSC204_getSensitivity(mK.T);
diff_oK = oK.sf - oK.sf_calc;
diff_mK = mK.sf - mK.sf_calc;

% plots Sensitivitäten
figNum = figNum + 1; figure(figNum); clf(figNum);
subplot(3,1,1); plot(t, oK.sf, '-', t, oK.sf_calc, '--');
title('(oK): Sensitivity factors in [uV/mW]'); legend('measured', 'calc''d'); timelabel(); ylabel('uV/mW')
subplot(3,1,2); plot(t, mK.sf, '-', t, mK.sf_calc, '--');
title('(mK): Sensitivity factors in [uV/mW]'); legend('measured', 'calc''d'); timelabel(); ylabel('uV/mW')
subplot(3,1,3); plot(t, diff_oK, t, diff_mK); title('Differences')

% Check: Sensitivitäten stimmen überein
sens_maxDiff_oK = max(abs(diff_oK));
sens_maxDiff_mK = max(abs(diff_mK));
fprintf('Sensitivity differences (maxabs):   oK: %g    mK: %g\n', sens_maxDiff_oK, sens_maxDiff_oK);
sens_tol = 0.0001;
if (sens_maxDiff_oK > sens_tol) || (sens_maxDiff_mK > sens_tol)
   warning('Sensitivity differing!');
end


%% POWER PLOTS in mW

% Berechne mW (benutze eigene Sensitivitätsfaktoren)
oK.mW = DSC204_uV_to_mW(oK.uV, oK.T);  
mK.mW = DSC204_uV_to_mW(mK.uV, mK.T);
fixT = 125; % feste Sensitivitäts-Temperatur in degC
oK.mW_fixT = DSC204_uV_to_mW(oK.uV, fixT);
mK.mW_fixT = DSC204_uV_to_mW(mK.uV, fixT);

% plot
figNum = figNum + 1; figure(figNum); clf(figNum);
subplot(2,1,1); plot(t, [oK.mW oK.mW_fixT], t, oK.mW-oK.mW_fixT);
title('(oK): Power in mW'); timelabel(); legend('measured T',sprintf('Tfix=%g',fixT), 'diff', 'Location', 'NorthWest');
grid on; axis tight
subplot(2,1,2); plot(t, [mK.mW mK.mW_fixT], t, mK.mW-mK.mW_fixT);
title('(mK): Power in mW'); timelabel(); legend('measured T',sprintf('Tfix=%g',fixT), 'diff', 'Location', 'NorthWest');
grid on; axis tight






%% export all the stuff to base workspace
varlist = who();
for k = 1:length(varlist)
   % break  % uncomment to skip
   var = varlist{k};
   try
      assignin('base', var, eval(var));
   end
end


%% finito
return











%% HELPERS
   function timelabel()
      xlabel('time in min');
   end

   function ax2 = newaxes()
      axh = gca;
      ax2 = axes('Position',get(axh,'Position'),'XAxisLocation','bottom','YAxisLocation','right','Color','none','nextplot','add');
   end

end




