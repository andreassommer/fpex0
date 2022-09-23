% Test for frasersuzuki initial condition 
% -> compares FD derivative approximation to the derivatives from transformed symbolic AD.

r  =   2;
h  =  40;
z  = 135.0;
wr =  15;
sr =   0.1;
dx = 0.001;
x  = 3:dx:250;
x  = 0:0.001:137;
x = reshape(x,[],1); % ensure x is a single column

p = [r; h; z; wr; sr];

fsfunction = @frasersuzuki;           % use non-symbolic version
%fsfunction = @frasersuzuki_symbolic;  % use symbolic version

% Funktionswerte
y = fsfunction(x,p);

% Symbolische Ableitung auswerten
[yy, dydp, dydx] = fsfunction(x,p); 
dydr  = dydp(:,1);
dydh  = dydp(:,2);
dydz  = dydp(:,3);
dydwr = dydp(:,4);
dydsr = dydp(:,5);

figure(6231)
clf 
disp('Showing analytical derivatives')
subplot(2,5, 1); hold on; plot(x, y    , '-'); ylabel('y')
subplot(2,5, 5); hold on; plot(x, dydx , '-'); ylabel('dydx')
subplot(2,5, 6); hold on; plot(x, dydr , '-'); ylabel('dydr')
subplot(2,5, 7); hold on; plot(x, dydh , '-'); ylabel('dydh')
subplot(2,5, 8); hold on; plot(x, dydz , '-'); ylabel('dydz')
subplot(2,5, 9); hold on; plot(x, dydsr, '-'); ylabel('dydsr')
subplot(2,5,10); hold on; plot(x, dydwr, '-'); ylabel('dydwr')
drawnow


% Finite difference approximation
FDp_h = 1.0d-8;
FDp_dirs = eye(5);
y2r  = fsfunction(x,p+FDp_h*FDp_dirs(:,1));
y2h  = fsfunction(x,p+FDp_h*FDp_dirs(:,2));
y2z  = fsfunction(x,p+FDp_h*FDp_dirs(:,3));
y2wr = fsfunction(x,p+FDp_h*FDp_dirs(:,4));
y2sr = fsfunction(x,p+FDp_h*FDp_dirs(:,5));
dydr_FD  = (y2r-y)/FDp_h;
dydh_FD  = (y2h-y)/FDp_h;
dydz_FD  = (y2z-y)/FDp_h;
dydwr_FD = (y2wr-y)/FDp_h;
dydsr_FD = (y2sr-y)/FDp_h;

FDx_h = 1e-6;
y2x  = fsfunction(x+FDx_h, p);
dydx_FD = (y2x-y)/FDx_h;



pause
clf;
disp('Showing analytical derivatives and FD derivative approximations')
subplot(2,5, 1); hold on; plot(x, y    , '-'); ylabel('y')
subplot(2,5, 5); hold on; plot(x, dydx , '-'); plot(x, dydx_FD , '.'); ylabel('dydx')
subplot(2,5, 6); hold on; plot(x, dydr , '-'); plot(x, dydr_FD , '.'); ylabel('dydr')
subplot(2,5, 7); hold on; plot(x, dydh , '-'); plot(x, dydh_FD , '.'); ylabel('dydh')
subplot(2,5, 8); hold on; plot(x, dydz , '-'); plot(x, dydz_FD , '.'); ylabel('dydz')
subplot(2,5, 9); hold on; plot(x, dydsr, '-'); plot(x, dydsr_FD, '.'); ylabel('dydsr')
subplot(2,5,10); hold on; plot(x, dydwr, '-'); plot(x, dydwr_FD, '.'); ylabel('dydwr')
drawnow


pause
clf
disp('Showing difference between analytical derivatives and FD approximations')
subplot(2,5, 5); hold on; plot(x, dydx-dydx_FD , '.-'); ylabel('dydx')
subplot(2,5, 6); hold on; plot(x, dydr-dydr_FD , '.-'); ylabel('dydr')
subplot(2,5, 7); hold on; plot(x, dydh-dydh_FD , '.-'); ylabel('dydh')
subplot(2,5, 8); hold on; plot(x, dydz-dydz_FD , '.-'); ylabel('dydz')
subplot(2,5, 9); hold on; plot(x,dydsr-dydsr_FD, '.-'); ylabel('dydsr')
subplot(2,5,10); hold on; plot(x,dydwr-dydwr_FD, '.-'); ylabel('dydwr')
drawnow
