% Simple tester for FokkerPlanckODE

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
u0 = normpdf(xgrid, mu, sigma);
p  = [];
driftFcn     = @(t,p) -0.00008 * t ;
diffusionFcn = @(t,p) 0.000001 * t;

% time grid
tgrid = 1:2:100;
t0tf = tgrid([1 end]);

% generate function and jacobian handle
rhsfun   = @(t,u) FokkerPlanckODE(t, u, p, h, driftFcn, diffusionFcn);
Jacobian = @(t,u) FokkerPlanckODE_dfdu(t, u, p, h, driftFcn, diffusionFcn);

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
             
