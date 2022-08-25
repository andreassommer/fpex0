N = 1001;     % number of grid points
h = 1/(N-1);  % interval length

xl = 90;
xr = 250;
xgrid = linspace(xl,xr,N);

mu    = 130;
sigma = 3.5;
u0 = normpdf(xgrid, mu, sigma);
p  = [];
driftFcn     = @(t,p) -0.003;
diffusionFcn = @(t,p) 0.00001;

tspan = 1:100;

opts = odeset('AbsTol', 1e-6, 'RelTol', 1e-10, 'BDF', 'on', 'stats', 'on', 'vectorized', 'on');

profile clear; profile on
tic
sol = ode15s(@(t,u) FokkerPlanck_ODE(t, u, p, h, driftFcn, diffusionFcn), tspan, u0, opts);
toc
profile off
profile viewer

u = deval(sol, tspan);

figure(337); clf; hold('on')
for i = 1:length(tspan)
   tvec = tspan(i) * ones(N,1);
   plot3(xgrid, tvec, u(:,i));
end
view([6,27])
ylabel('Time t')
xlabel('Space x')
zlabel('Value u')
             
