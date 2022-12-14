function FPEX0_visualize(FPEX0setup, sol)
% FPEX0_visualize(FPEX0setup, sol)
%
% Visualizes a simulation sol generated by FPEX0_simulate.
%
% INPUT: 
%     FPEX0setup --> FPEX0 setup object
%            sol --> solution object (from Matlab integrator)    
%
% OUTPUT   none
%
% Author:  Andreas Sommer, Aug-Sep2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu


% accessors
fignum = 337;
tgrid = FPEX0setup.Grid.gridTdot;
xgrid = FPEX0setup.Grid.gridT;
N = FPEX0setup.Grid.N;

% evaluation
u = deval(sol, tgrid);

% plot
figure(fignum);
clf; hold('on')
for i = 1:length(tgrid)
   tvec = tgrid(i) * ones(N,1);
   plot3(xgrid, tvec, u(:,i));
end
view([6,27])
ylabel('Heat rate \beta  (time t)')
xlabel('Temperature T   (space x)')
zlabel('c_p   (value u)')
grid('on')

% if measurements are available, add them
for m = FPEX0setup.Measurements
   plot3(m.temperatures, m.heatrate, m.values, 'r.-');
end

end