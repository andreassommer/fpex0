%function test_VDE()

% simple per for VDE

FPEX0setup = FPEX0_exampleSetup();
FPEX0setup.Integration.options.Stats = 'on';

gridT   = FPEX0setup.Grid.gridT;
t0tf    = FPEX0setup.Grid.gridTdot([1 end]);
betamax = FPEX0setup.betamax;

% comparison grid
compCount       = 10;
compGridIndices = round(linspace(gridT(1), gridT(end), compCount)); 
compGridTdot    = linspace(0, betamax, 100);


% fraser suzuki parametrization
r  =     2;
h  =    40;
z  = 135.0;
wr =    15;
sr =   0.1;
p_IC = [r h z wr sr];

% FP drift and diffusion parametrization
p_FPdrift     = [  1.0    -0.0001   ];
p_FPdiffusion = [  1.001   0.001   ];

% assemble p
pvec = [ p_FPdrift  p_FPdiffusion  p_IC];

% dimension stuff
np    = length(pvec);           % number of parameters
np_FP = length(p_FPdrift) + length(p_FPdiffusion);
np_IC = length(p_IC);
N     = FPEX0setup.Grid.N;      % number of states


%% nominal solution
disp('===== NOMINAL =====')
tic();
solNOMINAL = FPEX0_simulate(FPEX0setup, pvec);
tt = toc();
fprintf('--- Took %gs\n', tt);
u_NOMINAL = deval(solNOMINAL, compGridTdot);



%% prepare the VDE

% get VDE integrator and options
VDEintegrator = FPEX0setup.Integration.VDEintegrator;
VDEoptions    = FPEX0setup.Integration.VDEoptions;

h       = FPEX0setup.Grid.h;
betamax = FPEX0setup.betamax;

% initial condition for VDE
[u0, Gp_IC] = FPEX0setup.IniDistFcn(gridT, p_IC);  % initial distribution u0 with derivative du0/dp0
Gp_FP = zeros(N, np_FP);       % for FP-parameters, no dependency of u0 on pFP
Gp0   = [ Gp_FP , Gp_IC ];     % initial value GP(0)
Gp0   = reshape(Gp0, [], 1);   % reshape into vector

FP_VDE_rhs = @(t,x) FokkerPlanckVDE(t, solNOMINAL, x, h, p_FPdrift, p_FPdiffusion, betamax, np, false);
FP_VDE_jac = @(t,x) FokkerPlanckVDE(t, solNOMINAL, x, h, p_FPdrift, p_FPdiffusion, betamax, np, true); % NOT IMPLEMENTED

%%
disp('===== VDE =====')
tic
solVDE = VDEintegrator(FP_VDE_rhs, t0tf, Gp0, VDEoptions);
tt = toc();
fprintf('--- Took %gs\n', tt);

% evaluate VDE
uu_VDE = deval(solVDE, compGridTdot);
for i = 1:np
   psensVDE{i} = uu_VDE((i-1)*N+1:i*N, :);  %#ok<SAGROW>
end



%% nominal solution and VDE using simulate
disp('===== NOMINAL and VDE with simulate() =====')
tic();
[XsolNOMINAL, XsolSENS] = FPEX0_simulate(FPEX0setup, pvec);
tt = toc();
fprintf('--- Took %gs\n', tt);
% nominal solution should be identical
fprintf('nominal identity:  %d\n', isequal(solNOMINAL.y, XsolNOMINAL.y));
Xuu_VDE = deval(XsolSENS, compGridTdot);
for i = 1:np
   XpsensVDE{i} = Xuu_VDE((i-1)*N+1:i*N, :);  %#ok<SAGROW>
end



%% make FD approximation for parameter sensitivities
disp('===== FD =====')
for i = 1:np
   [pFD,hFD] = FDdisturb(pvec, i);
   tic();
   solFDp{i}  = FPEX0_simulate(FPEX0setup, pFD);  %#ok<SAGROW>
   u_FD{i}    = deval(solFDp{i}, compGridTdot);   %#ok<SAGROW>
   psensFD{i} = ( u_FD{i} - u_NOMINAL ) / hFD;    %#ok<SAGROW>
   tFD = toc();
   fprintf('FD %d/%d took %gs\n', i, np, tFD);
end



%% compare
figbase = 680;
for i = 1:np
   figure(figbase); clf; hold on;
   for k = 1:compCount
      subplot(compCount, 1, k); hold on;
      plot(psensFD{i}  (compGridIndices(k),:), 'r' , 'DisplayName', 'FD' );
      plot(psensVDE{i} (compGridIndices(k),:), 'b' , 'DisplayName', 'VDE');
      plot(XpsensVDE{i}(compGridIndices(k),:), 'c:', 'DisplayName', 'VDE_X');
      legend('show')
   end
   pause
end
% 
% disp(psens)
% 
% % display solution
% FPEX0_visualize(FPEX0setup, sol);






% finito
return


% Helper
function [pp, h] = FDdisturb(p,i)
   if abs(p(i)) < 1e-10
      h = 1e-6;
   else
      h = p(i)*1e-6;
   end
   p(i) = p(i) + h;
   pp = p;
end

%end % function