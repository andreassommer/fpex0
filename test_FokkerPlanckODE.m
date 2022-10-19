function varargout = test_FokkerPlanckODE()
% Simple tester for FokkerPlanckODE
%#ok<*NBRAK>  


% number of grid points and interval length
N = 2001;
h = 1/(N-1);

% Grid
xl = 90;
xr = 280;
xgrid = linspace(xl,xr,N);

% generate sparsity pattern of jacobian
e = ones(N,1);
JPattern = spdiags([e e e], -1:1, N, N); 

% initial distribution, drift, diffusion
mu    = 130;
sigma = 3.5;
u0    = reshape( normpdf(xgrid, mu, sigma) , [], 1 );
driftParams     = [ 1e-6  -1e-4 ]; 
diffusionParams = [ 1e-6  1e-4 ];

% time grid ( = heating rates)
betamax = 100;
tgrid = 1:2:betamax;
t0tf = tgrid([1 end]);

% generate function and jacobian handle
rhsfun   = @(t,u) FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax);
Jacobian = @(t,u) FokkerPlanckODE(t, u, h, driftParams, diffusionParams, betamax, false, true, false);

% test/compare jacobians
% JacobianOLD = @(t,u) old_FokkerPlanckODE_dfdu(t, u, h, [], driftParams, [], diffusionParams);
% snorm = @(j) norm(j(:));
% jOLD = JacobianOLD(0, u0);
% jNEW = Jacobian(0, u0);
% jDIFF = full(max(max(abs(jOLD-jNEW))));
% fprintf('\n||jOLD||=%g  ||jNEW||=%g   Difference in Jacobian: %g\n\n', snorm(jOLD), snorm(jNEW), jDIFF);


% options and integrator selection
opts = odeset( 'AbsTol'      , 1e-4     ...
             , 'RelTol'      , 1e-8     ...
             , 'BDF'         , 'off'    ...
             , 'vectorized'  , 'on'     ...
             , 'JPattern'    , JPattern ...
             , 'Jacobian'    , Jacobian ...
       ...   , 'NonNegative' , 1:N      ...
             , 'stats'       , 'on'      );
             
integrator = @ode15s;   % use an implicit solver! ode15s or ode23s

% integration and measurements
profile clear; profile on
tic
sol = integrator(rhsfun, t0tf, u0, opts);
toc
profile off
profile viewer

% evaluation
u = deval(sol, tgrid);

% plot
figure(337); clf; hold('on')
for i = 1:length(tgrid)
   tvec = tgrid(i) * ones(N,1);
   plot3(xgrid, tvec, u(:,i));
end
view([6,27])
ylabel('Time t')
xlabel('Space x')
zlabel('Value u')
        
% output requestes?
if (nargout > 0)
   varargout{1} = u;
end

% finito
return

end
