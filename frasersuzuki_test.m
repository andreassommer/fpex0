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


% Funktionswerte
y = frasersuzuki(x,p);

% Symbolische Ableitung auswerten
[yy, dydp] = frasersuzuki(x,p); 
dydr  = dydp(:,1);
dydh  = dydp(:,2);
dydz  = dydp(:,3);
dydwr = dydp(:,4);
dydsr = dydp(:,5);

figure(6231)
clf 
disp('Showing analytical derivatives')
subplot(3,2,1); plot(x,y); ylabel('y')
subplot(3,2,2); hold on; plot(x,dydr , '-'); ylabel('dydr')
subplot(3,2,3); hold on; plot(x,dydh , '-'); ylabel('dydh')
subplot(3,2,4); hold on; plot(x,dydz , '-'); ylabel('dydz')
subplot(3,2,5); hold on; plot(x,dydsr, '-'); ylabel('dydsr')
subplot(3,2,6); hold on; plot(x,dydwr, '-'); ylabel('dydwr')
drawnow


% % Finite Differenzen Approximation
FDh = 1.0d-8;
FDdirs = eye(5);
y2r  = frasersuzuki(x,p+FDh*FDdirs(:,1));
y2h  = frasersuzuki(x,p+FDh*FDdirs(:,2));
y2z  = frasersuzuki(x,p+FDh*FDdirs(:,3));
y2wr = frasersuzuki(x,p+FDh*FDdirs(:,4));
y2sr = frasersuzuki(x,p+FDh*FDdirs(:,5));

dydr_FD  = (y2r-y)/FDh;
dydh_FD  = (y2h-y)/FDh;
dydz_FD  = (y2z-y)/FDh;
dydwr_FD = (y2wr-y)/FDh;
dydsr_FD = (y2sr-y)/FDh;

pause
clf;
disp('Showing analytical derivatives and FD derivative approximations')
subplot(3,2,1); plot(x,y); ylabel('y')
subplot(3,2,2); hold on; plot(x,dydr , '-'); plot(x,dydr_FD , '.'); ylabel('dydr')
subplot(3,2,3); hold on; plot(x,dydh , '-'); plot(x,dydh_FD , '.'); ylabel('dydh')
subplot(3,2,4); hold on; plot(x,dydz , '-'); plot(x,dydz_FD , '.'); ylabel('dydz')
subplot(3,2,5); hold on; plot(x,dydsr, '-'); plot(x,dydsr_FD, '.'); ylabel('dydsr')
subplot(3,2,6); hold on; plot(x,dydwr, '-'); plot(x,dydwr_FD, '.'); ylabel('dydwr')
drawnow


pause
clf
disp('Showing difference between analytical derivatives and FD approximations')
subplot(3,2,2); hold on; plot(x,dydr-dydr_FD , '.-'); ylabel('dydr')
subplot(3,2,3); hold on; plot(x,dydh-dydh_FD , '.-'); ylabel('dydh')
subplot(3,2,4); hold on; plot(x,dydz-dydz_FD , '.-'); ylabel('dydz')
subplot(3,2,5); hold on; plot(x,dydsr-dydsr_FD, '.-'); ylabel('dydsr')
subplot(3,2,6); hold on; plot(x,dydwr-dydwr_FD, '.-'); ylabel('dydwr')
drawnow
