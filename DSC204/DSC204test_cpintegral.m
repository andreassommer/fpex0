function [dsc, dsc_reference] = DSC204test_cpintegral(dsc, dsc_reference)
% JUST TESTING!


% if no dscdata is specified, load it
if (nargin < 2)
   dsc_reference = DSC204_readFiles('Sap-Kurve_10Kmin_H_Segment_7.csv');  % load saphire
end
if (nargin < 1)
   dsc = DSC204_readFiles('ExpDat_16-*_mitKorr_*Kmin_H.csv');         % load dsc measurements
end

% always re-calculate cp
dsc = DSC204_addCP(dsc, dsc_reference);                               % add cp information
dsc = DSC204_sortByTstep(dsc);                                        % sort by Tstep
dsc = DSC204_sortByDescField(dsc, 'IDENTITY');                        % sort by Indentity


% save data in main workspace
assignin('base','xxxdsc',dsc);
assignin('base','xxxdscref', dsc_reference);
   
% display and plot
for k = 1:length(dsc)
   
   % accessors to data
   data = dsc(k);
   T  = data.cp.T;
   cp = data.cp.fun(T);
   cpfun = data.cp.fun;
   bl = data.cp.blfun(T);
   cplatent = data.cp.latentfun(T);
   cplatentfun = data.cp.latentfun;
   
   % heat of fusion (integral over latentcp)
   deltaH  = integral(cplatentfun, T(1), T(end));
   
   % upper limit of heat of fusion (integral with linear baseline starting from mean(cp(T(1:w))) to mean(cp(T(end-w:end)));
   % window length: one degree celsius
   window = find(T-T(1)>=1.0);
   mean_cp_l = mean(cp(1:window));
   mean_cp_r = mean(cp(end-window:end));
   linbase = @(x) interp1([T(1), T(end)], [mean_cp_l, mean_cp_r], x);
   deltaH_ub = integral(@(x) max(0, cpfun(x)-linbase(x)), T(1), T(end));
   
   fprintf('Displaying #%02d: %-40s --- deltaH = %5.1f -- deltaH_ub = %5.1f -- Tonset = %5.1f -- Toffset = %5.1f -- Tmin/Tmax=(%5.1f,%5.1f)\n', ...
            k, data.fileSpec, deltaH, deltaH_ub, data.cp.bldata.onset, data.cp.bldata.endset, T(1), T(end));
   
   if true
      clf; hold on;
      plot(T,cp,'b');
      plot(T,bl,'g','linewidth',2.0);
      plot(T,cplatent,'r');
      plot(T,linbase(T),'g:');
      legend('cp', 'baseline', 'cp_{latent}','simple baseline','Location','West')
      titlestr = sprintf('File: %s --- Identity: %s', data.fileSpec, data.desc.IDENTITY);
      title(titlestr, 'Interpreter', 'none')
      ax = axis(); 
      infostr = makeInfo('deltaH'            , deltaH                       , ...             
                         'deltaH (simpleBL)' , deltaH_ub                    , ...
                         'Tonset'            , data.cp.bldata.onset         , ...
                         'Tendset'           , data.cp.bldata.endset        , ...
                         'Tmin'              , T(1)                         , ...
                         'Tmax'              , T(end)                       , ...
                         'left BL intercept' , data.cp.bldata.reg_l.a       , ...
                         'left BL slope'     , data.cp.bldata.reg_l.b       , ...
                         'left BL stddev'    , sqrt(data.cp.bldata.reg_l.s2), ...
                         'right BL intercept', data.cp.bldata.reg_r.a       , ...
                         'right BL slope'    , data.cp.bldata.reg_r.b       , ...
                         'right BL stddev'   , sqrt(data.cp.bldata.reg_r.s2) );
      text(ax(1), ax(4), infostr, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Interpreter', 'none', 'FontName', 'FixedWidth')
      pause
   end
   
end


   function infostr = makeInfo(varargin)
      maxlen = max(cellfun(@length, {varargin{1:2:end}}));
      formatstr = sprintf('\n %%%ds = %%7.3f', maxlen);
      infostr = '';
      for j=1:2:length(varargin)
         infostr = [infostr sprintf(formatstr, varargin{j}, varargin{j+1})];
      end
   end

% finito
end