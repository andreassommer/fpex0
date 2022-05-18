function datavec = DSC204_clearZeroFromMax(datavec)
   % datavec = DSC204_clearZeroFromMax(datavec)
   % 
   % Locates the maximum value in datavec, searches the first occurances
   % of a zero value to the left and to the right from the peak position
   % and deletes the noise before and after these positions.
   %
   % INPUT:  datavec --> data vector to be cleared
   %
   % OUTPUT: datavec --> cleared data vector
   %
   % Author:  Andreas Sommer, Apr2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %

   % Find index/position of maximal value
   [~, maxidx] = max(datavec);
   
   % Find the positions of the first zeros, seen from the peak value
   clearToIdx   = find(datavec(1:maxidx)==0, 1, 'last');
   clearFromIdx = find(datavec(maxidx:end)==0, 1, 'first') + maxidx;
   
   % delete everything that is farer away
   if ~isempty(clearToIdx),   datavec(1:clearToIdx)     = 0; end
   if ~isempty(clearFromIdx), datavec(clearFromIdx:end) = 0; end
end