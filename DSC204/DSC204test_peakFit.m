function [peakX, peakVal, peakRate] = DSC204test_peakFit(dscdata)

% config
fignum = 532;
linewidth = 2.0;

% sord by heat rate and group by ID
dscdata = DSC204_sortByTstep(dscdata, 'ascend');
dscgroups = DSC204_groupByID(dscdata);

% determine peaks
peakX = cell(size(dscgroups));
peakRate = cell(size(dscgroups));
peakVal = cell(size(dscgroups));

for grp = 1:length(dscgroups)
   
   for k = 1:length(dscgroups{grp})
      dsc = dscgroups{grp}(k);
      %X = dsc.data(:,1);  Y = dsc.data(:,3);      % DSC µV signals
      X = dsc.cp.T     ;  val = dsc.cp.latentdata;  % latent CP
      % determine peaks
      [maxVal, maxIdx] = max(val);
      maxX = X(maxIdx);
      % store peak info
      peakX{grp}(end+1) = maxX;
      peakVal{grp}(end+1) = maxVal;
      peakRate{grp}(end+1) = dsc.rate;
   end
   
   % plot
   figure(fignum);
   subplot(1, length(dscgroups), grp)
   plot(peakX{grp}, peakVal{grp}, '.--', 'LineWidth', linewidth)
   
   % try to fit polynomial
   
end


end