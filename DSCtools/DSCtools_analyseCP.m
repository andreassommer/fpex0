function DSCtools_analyseCP(DSCsample, varargin)
% DSCtools_analyseCP(DSCsample, varargin)
%
% Calculates latend heat of fusion (deltaH) with all available baselines and display results graphically.
%
% INPUT:    DSCsample --> dsc data structure with added cp information
%     key-value-pairs --> 'fignum'   for plotting window
%                         'fontsize' for setting the font size
% 
% OUTPUT:     none
%
%
% Author: Andreas Sommer, Jul2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu


% arg checks
if (nargin < 1), error('Need input.'); end

% process args
args = varargin;
[fignum   , args] = olGetOption(args, 'fignum'  , 913);
[fontsize , args] = olGetOption(args, 'fontsize',  []);
olWarnIfNotEmpty(args);

% check if we have a cp field
if ~isfield(DSCsample, 'cp')
   error('DSC data does not have a cp field. Use DSCtools_addCP() first.')
end

% retrieve information about latent heat of fusion
[deltaH, deltaHinfo] = DSCtools_calc_deltaH(DSCsample, 'doDisplay', true);

% display and plot
for k = 1:length(DSCsample)
   
   % accessors to data
   dk = DSCsample(k);
   T        = dk.cp.T;
   cp       = dk.cp.fun(T);
   cplatent = dk.cp.latentfun(T);
   bldata   = dk.cp.bldata;
   Tonset   = bldata.onset;
   Tendset  = bldata.endset;
   dHk      = deltaHinfo(k);  % accessor for deltaHinfo

   % short info
   makeMessage('#%02d: %-40s  dH = %5.2f    dH_ub = %5.2f   dH_ToTe = %5.2f    Tonset = %5.2f   Tendset = %5.2f   (Tmin,Tmax)=(%5.1f,%5.1f)\n', ...
                   k,  dk.ID, dHk.bl_linrange.deltaH, dHk.bl_simple.deltaH  ,dHk.bl_ToTe.deltaH, Tonset, Tendset, T(1), T(end));
   
   % if a fignum is specified, do a visualization
   if ~isempty(fignum)
      figure(fignum); clf; axh = gca(); hold(axh, 'on');
      plot(axh, T, cp                                   , 'c' , 'DisplayName', 'cp'           , 'linewidth', 1.5);
      plot(axh, T, cplatent                             , 'm' , 'DisplayName', 'cp latent'    , 'linewidth', 1.5);
      plot(axh, T, bldata.baselines.bl_linrange.blfun(T), 'g' , 'DisplayName', 'baseline'     , 'linewidth', 1.0);
      plot(axh, T, bldata.baselines.bl_simple  .blfun(T), 'g:', 'DisplayName', 'baseline_ub'  , 'linewidth', 1.0);
      plot(axh, T, bldata.baselines.bl_ToTe   .blfun(T), 'r:', 'DisplayName', 'baseline_ToTe', 'linewidth', 1.5);
      xline(axh, [bldata.onset  bldata.endset], 'r:' , 'onset/endset' , 'LabelVerticalAlignment', 'bottom', 'DisplayName', 'onset/endset');
      xline(axh, [bldata.linL   bldata.linR  ], 'r--', 'linear ranges', 'LabelVerticalAlignment', 'bottom', 'DisplayName', 'linear range');
      legend(axh, 'Location','NorthEast','Interpreter','none');
      titlestr = sprintf('#%d/%d  Identity: %s', k, length(DSCsample), dk.ID);
      title(axh, titlestr, 'Interpreter', 'none')
      infostr = makeInfo( 'deltaH'            , dHk.bl_linrange.deltaH  ...             
                        , 'deltaH (Ton-Tend)' , dHk.bl_ToTe.deltaH      ...
                        , 'deltaH (simpleBL)' , dHk.bl_simple.deltaH    ...
                        , 'Tonset'            , bldata.onset          ...
                        , 'Tendset'           , bldata.endset         ...
                        , 'linear limit left' , bldata.linL           ...
                        , 'linear limit right', bldata.linR           ...
                        , 'Tmin'              , T(1)                  ...
                        , 'Tmax'              , T(end)                ...
                        , 'BL left icept'     , bldata.reg_L.a        ...
                        , 'BL left slope'     , bldata.reg_L.b        ...
                        , 'BL left stddev'    , sqrt(bldata.reg_L.s2) ...
                        , 'BL right icept'    , bldata.reg_R.a        ...
                        , 'BL right slope'    , bldata.reg_R.b        ...
                        , 'BL right stddev'   , sqrt(bldata.reg_R.s2) ...
                        );
      axlim = axis(); 
      th = text(axlim(1), axlim(4), infostr, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Interpreter', 'none', 'FontName', 'FixedWidth');
      if ~isempty(fontsize)
         set([th axh], 'FontSize', fontsize); 
      end
      pause
   end
   
end % of visualization

% finito
return

end % of function


%% HELPERS
function infostr = makeInfo(varargin)
   maxlen = max(cellfun(@length, varargin(1:2:end)));
   formatstr = sprintf('\n %%-%ds = %%.4g', maxlen);
   infostr = '';
   for j=1:2:length(varargin)
      infostr = [infostr sprintf(formatstr, varargin{j}, varargin{j+1})];  %#ok<AGROW>
   end
end


