% number of grid points and interval length
N = 2001;
h = 1/(N-1);

% Grid
xl = 90;
xr = 250;
xgrid = linspace(xl,xr,N);

% generate sparsity pattern of jacobian
e = ones(N,1);
JPattern = spdiags([e e e], -1:1, N, N); 

% initial distribution, drift, diffusion
mu    = 130;
sigma = 3.5;
u0 = normpdf(xgrid, mu, sigma);
p  = [];
driftFcn     = @(t,p) -0.003;
diffusionFcn = @(t,p) 0.00001;

% time grid
tspan = 1:100;

% options and integrator selection
opts = odeset('AbsTol', 1e-6, 'RelTol', 1e-10, 'BDF', 'off', 'stats', 'on', 'vectorized', 'on', 'JPattern', JPattern);
integrator = @ode15s;

% integration and measurements
profile clear; profile on
tic
sol = integrator(@(t,u) FokkerPlanckODE(t, u, p, h, driftFcn, diffusionFcn), tspan, u0, opts);
toc
profile off
profile viewer

% evaluation
u = deval(sol, tspan);

% plot
figure(337); clf; hold('on')
for i = 1:length(tspan)
   tvec = tspan(i) * ones(N,1);
   plot3(xgrid, tvec, u(:,i));
end
view([6,27])
ylabel('Time t')
xlabel('Space x')
zlabel('Value u')
             
