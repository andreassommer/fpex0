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
%                            'group'      : one of 'heatrate' or 'mass'            [default: heatrate         ]
%                            'ngroups'    : expected number of groups              [default: [] for autodetect]
%                            'linkaxes'   : forwarded to matlab's linkaxes()       [default: 'off'            ]
%                            'grid'       : forwarded to matlab's grid()           [default: 'on'             ]
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
[normbymass  , args] = olGetOption(args, 'normbymass'  , true      );
[normbyrate  , args] = olGetOption(args, 'normbyrate'  , true      );
[linkopt     , args] = olGetOption(args, 'linkaxes'    , 'off'     );
[gridopt     , args] = olGetOption(args, 'grid'        , 'off'     );
olWarnIfNotEmpty(args);

% data sort - done now by grouping
% DSCdata = DSCtools_sortByMass(DSCdata, 'ascend');
% DSCdata = DSCtools_sortByHeatrate(DSCdata, 'ascend');

% prepare grouping
switch lower(group)
   case {'heatrate', 'heatrates', 'rate', 'rates'}
      groupValues = [DSCdata.rate];
   case {'mass', 'masses'}
      groupValues = [DSCdata.mass];
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
[idx, C] = kmeans(reshape(groupValues, [], 1), groupCount, 'OnlinePhase', 'on', 'Start', 'uniform', 'replicates', 100);
[C, ord] = sort(C, 'ascend');
groupMeans = C;
groupIdx = zeros(size(idx));
for i = 1:numel(ord)
    groupIdx(idx == ord(i)) = i;
end

% prepare figure
handles.axes   = zeros(length(groupCount), 1);
handles.figure = figure(fignum);
clf(handles.figure, 'reset');
handles.layout = tiledlayout(handles.figure,'flow','TileSpacing','compact');

% prepare storage for plot handles
nDSC = length(DSCdata);
handles.plots = zeros(nDSC, 1);
plotidx = 0;

% plot by group
for k = 1:groupCount
   plotidx = plotidx + 1;           % for enumerating the plots in same order as DSCdata
   handles.axes(k) = nexttile(); 
   hold(handles.axes(k), 'on');
   sel = DSCdata(groupIdx == k);    % select those that belong to the k-th group
   for i = 1:length(sel)
      d = sel(i);
      T = d.data.T;
      Y = d.data.(datafield);
      if (normbymass),  Y = Y ./ d.mass;  end
      if (normbyrate),  Y = Y ./ d.rate;  end
      handles.plots(plotidx) = plot(handles.axes(k), T, Y, 'Displayname', makeTitleString(d, 'compact'));
   end
   legend(handles.axes(k), 'Interpreter', 'none', 'location', 'southeast')
   title(handles.axes(k), sprintf('Group #%d:  %s about %g', k, group, groupMeans(k)), 'interpreter', 'none');
end
unifyAxes(handles.axes, [], [], [], linkopt, gridopt);

% finito
if (nargout >= 1), varargout{1} = handles; end
return

end % of function




%% helper for title string generation

function str = addString(str, add)
   str = horzcat(str, '   ', add);
end

function tstring = makeTitleString(dk, type)
   file = dk.rawData.fileSpec;
   rate = dk.rate;
   mass = dk.mass;
   rep = str2double(extractAfter(dk.ID,'rep'));
   switch lower(type)
      case 'full'
         tstring = '';
         tstring = addString(tstring, sprintf('rate: %+g K/min', rate));
         tstring = addString(tstring, sprintf('mass: %g µg', mass));
         tstring = addString(tstring, sprintf('rep #%d', rep));
         tstring = addString(tstring, sprintf('filename: %s', file));
      case 'compact'
         tstring = sprintf('%+g K  %g µg  #%d', rate, mass, rep);
   end
end

