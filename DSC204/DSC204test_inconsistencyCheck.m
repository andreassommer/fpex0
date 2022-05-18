function DSC204test_inconsistencyCheck(oK, mK)

  % Keine Daten? Dann lade Daten ohne und mit Korrektur
  % Bedingung: Alle DSC-Messdaten liegen zusammen im aktuellen Verzeichnis!
  %   oK = DSC204_loadall('*ohneKorr*_H.csv');
  %   mK = DSC204_loadall('*mitKorr*_H.csv');
  
  warning('Handle this tool with care. Just for internal testing!')
  
  if (nargin<2)
     disp('Lade Datensatz.')
     disp('Bedingung: Alle DSC-Messdaten liegen zusammen im aktuellen Verzeichnis!')
     oK = DSC204_loadall('*ohneKorr*_H.csv');
     mK = DSC204_loadall('*mitKorr*_H.csv');
     fprintf('\nExportiere oK und mK in base-workspace\n')
     assignin('base', 'oK', oK);
     assignin('base', 'mK', mK);
     fprintf('\nSchneller Aufruf ab jetzt mit: DSC204_inconsistencyCheck(oK, mK)\n')
  end
  
  % choice
  selection = 1:6*7;  % 7 verschiedene Temperaturgradienten, 6 Experimente 407,408,409,416,417,418 
  
  
  % Zeige die Korrekturen an
  Tfig = figure(414); clf(Tfig); axT = gca; hold on;
  tfig = figure(415); clf(tfig); axt = gca; hold on; 
  for n = selection
     % Weil einige Dateien nicht das enthalten, was sie vom Dateinamen her versprechen,
     % muessen wir das in einen Try-Block wickeln um diese zu ueberspringen
     try 
        T = oK(n).data(:,1);                      % Temperatur
        t = oK(n).data(:,2);                      % Zeit
        diff = oK(n).data(:,3) - mK(n).data(:,3); % Korrektur = Differenz zwischen oK und mK
        plot(axT,T,diff);                         % Plot mit x-Achse = Temperatur
        plot(axt,t,diff);                         % Plot mit x-Achse = Zeit
        label = sprintf('#%d: %s', n, mK(n).fileSpec);
        if isfield(mK(n).desc, 'CORRFILE')
           label = sprintf('%s with CORRFILE: %s', label, mK(n).desc.CORRFILE);
        end
        idx = floor(length(T)/4);
        hT = text(T(idx), diff(idx), label, 'Parent', axT, 'Interpreter', 'none', 'Fontsize', 8);
        ht = text(t(idx), diff(idx), label, 'Parent', axt, 'Interpreter', 'none', 'Fontsize', 8);
     catch err
        fprintf('Ignoring error @ #%d for file %s: %s\n', n, mK(n).fileSpec, err.message)
     end
  end
  title(axT,'DSC-Korrektur ueber Temperatur aufgetragen.'); xlabel(axT,'Temperatur in degC'); ylabel(axT,'Korrektur');
  title(axt,'DSC-Korrektur ueber Zeit aufgetragen.')      ; xlabel(axt,'Zeit in min');        ylabel(axt,'Korrektur');
  
  % finito
end