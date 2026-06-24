function varargout = DSCtools_display(DSCdata, varargin)
% handles = DSCtools_display(DSCdata, key-value-pairs*)
%
% Displays the DSCdata, sorted by heat rate and mass, in individual plots.
%
% INPUT:         DSCdata --> DSC data structure as returned by DSCtools_readFile
%        key-value-pairs --> possible keys are:
%                            'fignum'     : figure to display into                 [default: 1000             ]
%                            'datafield'  : name of data field to be displayed     [default: mW               ]
%                            'normbymass' : if true, data is divided by mass       [default: true             ]
%                            'normbyrate' : if true, data is divided by heat rate  [default: true             ]
%                            'sort'       : one of 'heatrate' or 'mass' or 'none'  [default: 'none'           ]
%                            'group'      : one of 'heatrate' or 'mass' or 'none'  [default: heatrate         ]
%                            'ngroups'    : expected number of groups              [default: [] for autodetect]
%                            'linkaxes'   : forwarded to matlab's linkaxes()       [default: 'off'            ]
%                            'grid'       : forwarded to matlab's grid()           [default: 'on'             ]
%                            'multifig'   : plot groups in individual figures      [default: false            ]
%                            'clearfig'   : clear figure before plotting           [default: true             ]
%                            'mergefig'   : add to existing figure/layout          [default: true             ]
%
% OUTPUT:  handles --> structure with handles of created plot
%
% NOTE:  if the ngroups is not delivered, an automated clustering is performed:
%        - by the number of heat rates if 'group' is set to 'heatrate'
%        - by the number of masses if 'group' is set to 'mass'
%
% Author: Andreas Sommer, Jun2026
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu


% quick return on empty data
if isempty(DSCdata)
   warnID  = 'DSCTools:DSC_NO_DATA';
   warnMSG = 'Empty DSCdata structure';
   warning(warnID, warnMSG);
   return;
end

% args and defaults
args = varargin;
[fignum      , args] = olGetOption(args, 'fignum'      , 1000      );
[datafield   , args] = olGetOption(args, 'datafield'   , 'mW'      );
[group       , args] = olGetOption(args, 'group'       , 'heatrate');
[ngroups     , args] = olGetOption(args, 'ngroups'     , []        );
[sorting     , args] = olGetOption(args, 'sort'        , 'none'    );
[normbymass  , args] = olGetOption(args, 'normbymass'  , true      );
[normbyrate  , args] = olGetOption(args, 'normbyrate'  , true      );
[linkopt     , args] = olGetOption(args, 'linkaxes'    , 'off'     );
[gridopt     , args] = olGetOption(args, 'grid'        , 'off'     );
[multifig    , args] = olGetOption(args, 'multifig'    , false     );
[clearfig    , args] = olGetOption(args, 'clearfig'    , true      );
[mergefig    , args] = olGetOption(args, 'mergefig'    , true      );
[colors      , args] = olGetOption(args, 'colors'      , []        );
olWarnIfNotEmpty(args);

% inits
nDSC = length(DSCdata);

% if colors was given as scalar, replicate it
if istext(colors) && isscalar(colors)           ... % case 1: scalar string or char array like 'b' for blue
   || isnumeric(colors) && length(colors)==3        % case 2: single color rgb spec like [0.3 0.3 0.3] for dark gray
   colors = repmat(colors, nDSC, 1);
end

% sort if user requested - must be before grouping
switch lower(sorting)
   case {'heatrate', 'heatrates', 'rate', 'rates'}
      DSCdata = DSCtools_sortByHeatrate(DSCdata, 'ascend');
   case {'mass', 'masses'}
      DSCdata = DSCtools_sortByMass(DSCdata, 'ascend');
   case {'none'}
      % do nothing
   otherwise
      error('Unknown sorting request: %s', sorting)
end

% prepare grouping
switch lower(group)
   case {'heatrate', 'heatrates', 'rate', 'rates'}
      groupValues = [DSCdata.rate];
   case {'mass', 'masses'}
      groupValues = [DSCdata.mass];
   case {'none'}
      groupValues = [];
      % do nothing
   otherwise
      error('Unknown group: %s', group);
end

% use user specified group count
if isempty(ngroups)
   groupCount = length(unique(groupValues));
else
   groupCount = ngroups;
end

% assign a group to every DSCsample, sort by group value
if isempty(groupValues)
   groupCount = 1;
   groupMeans = nan();
   groupIdx = ones(nDSC, 1);
else
   [idx, C] = kmeans(reshape(groupValues, [], 1), groupCount, 'OnlinePhase', 'on', 'Start', 'uniform', 'replicates', 100);
   [C, ord] = sort(C, 'ascend');   % Matlab always uses stable sort
   groupMeans = C;
   groupIdx = zeros(size(idx));
   for i = 1:numel(ord)
      groupIdx(idx == ord(i)) = i;
   end
end

% prepare figure(s)
handles.axes   = zeros(length(groupCount), 1);
if multifig
   handles.figure = zeros(length(groupCount), 1);
   handles.layout = [];
else
   handles.figure = prepareFigure(fignum, clearfig);
   handles.layout = findall(handles.figure, 'Type', 'TiledLayout');
   if isempty(handles.layout)
      handles.layout = tiledlayout(handles.figure,'flow','TileSpacing','compact');
   end
end

% prepare storage for plot handles
handles.plots = zeros(nDSC, 1);
plotidx = 0;

% legend font settings
legendstyles = {'FontName', 'fixedwidth', 'FontSize', 8};

% use axis offset if mergefig is NOT active
axisoffset = 0;
if ~mergefig
   axisoffset = length(findall(handles.layout, 'type', 'axes'));
end

% plot by group
for k = 1:groupCount
   if multifig
      handles.figures(k) = prepareFigure(fignum+k, clearfig);      % in multifig, we start at +1 to avoid overwriting figure fignum
      handles.axes(k) = gca();
   else
      handles.axes(k) = nexttile(handles.layout, k+axisoffset);    % merge into k-th axis
   end
   hold(handles.axes(k), 'on');      % hold the axis for multiple plots
   sel = DSCdata(groupIdx == k);     % select those that belong to the k-th group
   for i = 1:length(sel)
      plotidx = plotidx + 1;         % for enumerating the plots in same order as DSCdata
      d = sel(i);
      T = d.data.T;
      Y = d.data.(datafield);
      if (normbymass),  Y = Y ./ d.mass;  end
      if (normbyrate),  Y = Y ./ d.rate;  end
      ph = plot(handles.axes(k), T, Y, 'Displayname', makeLabel(d, 'compact'));  % plot
      if ~isempty(colors)
         set(ph, 'color', colors(plotidx,:));
      end
      addDataTip(ph, d)
      handles.plots(plotidx) = ph;   % store the plot handle
   end
   legend(handles.axes(k), legendstyles{:}, 'Interpreter', 'none', 'location', 'southeast')
   title(handles.axes(k), sprintf('Group #%d:  %s about %g', k, group, groupMeans(k)), 'interpreter', 'none');
   % add tooltip

end
unifyAxes(handles.axes, [], [-inf inf], [], linkopt, gridopt);

% finito
if (nargout >= 1), varargout{1} = handles; end
return

end % of function




%% HELPERS

function str = addString(str, add)
   str = horzcat(str, '   ', add);
end

function tstring = makeLabel(dk, type)
   file = dk.rawData.fileSpec;
   rate = dk.rate;
   mass = dk.mass;
   switch lower(type)
      case 'full'
         tstring = '';
         tstring = addString(tstring, sprintf('rate: %+4g K/min', rate));
         tstring = addString(tstring, sprintf('mass: %5g mg', mass));
         tstring = addString(tstring, sprintf('filename: %s', file));
      case 'compact'
         tstring = sprintf('%+4g K  %5.2f mg  %s', rate, mass, file);
   end
end

function fh = prepareFigure(fignum, clearfig)
   fh = figure(fignum);
   if clearfig
      clf(fh, 'reset');
   end
end

function addDataTip(ph, dscdata)
   len = length(ph.XData);
   repstr = @(s) repmat(string(s), len, 1);
   repvec = @(v) repmat(v, len, 1);
   ph.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Mass'  , repvec(dscdata.mass));
   ph.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Rate'  , repvec(dscdata.rate));
   ph.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('ID'  , repstr(dscdata.ID));
   ph.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('File', repstr(dscdata.rawData.fileSpec));
end
