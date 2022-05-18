function blds = DSC204_getBaselineDetectionSettings()
% blds = DSC204_getBaselineDetectionSettings()
%
% Returns defaults for detection of linear ranges.
%
% INPUT:   none 
%          (this might be extended to choose from a preset choice)
%
% OUTPUT:  bldd --> structure containing BaseLine Detection Settings
%
%
% Author:  Andreas Sommer, Apr2017, Dec2021
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%

 % LEFT SIDE:               % absolute and relative deviations
 blds.L.absdevA  = -1;      % abs dev from initial y-axis intercept
 blds.L.absdevB  = -1;      % abs dev from initial slope
 blds.L.absdevS2 = -1;      % abs dev from initial standard deviation
 blds.L.reldevA  = inf;     % rel dev from initial y-axis intercept
 blds.L.reldevB  = 0.01;    % rel dev from initial slope
 blds.L.reldevS2 = 2.00;    % rel dev from initial standard deviation
 blds.L.initfraction = 0.1; % part/percentage for initial estimate

 % RIGHT SIDE:              % absolute and relative deviations
 blds.R.absdevA  = -1;      % abs dev from initial y-axis intercept
 blds.R.absdevB  = 0.05;    % abs dev from initial slope
 blds.R.absdevS2 = -1;      % abs dev from initial standard deviation
 blds.R.reldevA  = inf;     % rel dev from initial y-axis intercept
 blds.R.reldevB  = inf;     % rel dev from initial slope
 blds.R.reldevS2 = 2.00;    % rel dev from initial standard deviation
 blds.R.initfraction = 0.2; % part/percentage for initial estimate
 
 

 % SET FOR FPEX0 initial submission
%  tmpblds.absdevA  = -1;      % abs dev from initial y-axis intercept
%  tmpblds.absdevB  = 0.01;    % abs dev from initial slope
%  tmpblds.absdevS2 = -1;      % abs dev from initial standard deviation
%  tmpblds.reldevA  = inf;     % rel dev from initial y-axis intercept
%  tmpblds.reldevB  = inf;     % rel dev from initial slope
%  tmpblds.reldevS2 = 2.00;    % rel dev from initial standard deviation
%  tmpblds.initfraction = 0.2; % part/percentage for initial estimate
%  blds.R = tmpblds;
%  blds.L = tmpblds;
 
 
end