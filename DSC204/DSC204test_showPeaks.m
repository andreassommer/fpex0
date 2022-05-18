function DSC204test_showPeaks(dscdata)

fignum=123;
figure(fignum);
clf; 
subplot(2,1,1); hold on
subplot(2,1,2); hold on
groupeddsc = DSC204_groupByID(dscdata);

n = length(groupeddsc);

T     = cell(n,1);
cpmax = cell(n,1);
rates = cell(n,1);

styles = {'b.-', 'm.-', 'k.-', 'g.-', 'r.-', 'c.-'};

for k = 1:length(groupeddsc)
   dsc = groupeddsc{k};
   for j = 1:length(dsc)
      [maxval, maxidx] = max(dsc(j).cp.values);
      cpmax{k}(j) = maxval;
      T{k}(j) = dsc(j).cp.T(maxidx);
      rates{k}(j) = dsc(j).Tinfo.Tstep;
   end
   
   subplot(2,1,1);
   plot(T{k}, cpmax{k}, styles{mod(k,length(styles))}, 'DisplayName', dsc(1).desc.IDENTITY)
   
   subplot(2,1,2);
   plot(rates{k}, T{k}, styles{mod(k,length(styles))}, 'DisplayName', dsc(1).desc.IDENTITY)
end

subplot(2,1,1); xlabel('Temperature T'); ylabel('maximum cp value'); legend show
subplot(2,1,2); xlabel('rate dT/dt'); ylabel('Temperature at max cp'); legend show



end