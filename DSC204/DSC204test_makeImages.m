
% some setup
colors = {'b', 'g', 'r', 'c', 'm', 'y', 'k'};

% lade daten experimentweise
filepath = DSC204_getDSCDataPath();
disp('Lade Daten');
DSC407 = DSC204_readFiles(fullfile(filepath, 'ExpDat_16-407-3_mitKorr_*_H.csv'));
DSC408 = DSC204_readFiles(fullfile(filepath, 'ExpDat_16-408-3_mitKorr_*_H.csv')); 
DSC409 = DSC204_readFiles(fullfile(filepath, 'ExpDat_16-409-3_mitKorr_*_H.csv')); 
DSC416 = DSC204_readFiles(fullfile(filepath, 'ExpDat_16-416-3_mitKorr_*_H.csv')); 
DSC417 = DSC204_readFiles(fullfile(filepath, 'ExpDat_16-417-3_mitKorr_*_H.csv')); 
DSC418 = DSC204_readFiles(fullfile(filepath, 'ExpDat_16-418-3_mitKorr_*_H.csv')); 

% sortiere
if true
   DSC407 = DSC204_sortByTstep(DSC407);
   DSC408 = DSC204_sortByTstep(DSC408);
   DSC409 = DSC204_sortByTstep(DSC409);
   DSC416 = DSC204_sortByTstep(DSC416);
   DSC417 = DSC204_sortByTstep(DSC417);
   DSC418 = DSC204_sortByTstep(DSC418);
end

% prozessiere daten
if true
   DSC407 = DSC204_testprocessing(DSC407);
   DSC408 = DSC204_testprocessing(DSC408);
   DSC409 = DSC204_testprocessing(DSC409);
   DSC416 = DSC204_testprocessing(DSC416);
   DSC417 = DSC204_testprocessing(DSC417);
   DSC418 = DSC204_testprocessing(DSC418);
end


% plot
figure(998); clf;
subplot(2,3,1); plot(DSC407.T, DSC407.uVn)


