% parameters
% a0 --> peak area
% a1 --> center of Gaussian component
% a2 --> standard deviation of the Gaussian component
% a3 --> peak distortion
p = [198 ; 120 ; 5 ; -1];
p = [228  ;  110  ;  3 ; -2.5];

% x range
x = 60:0.1:150;
x = reshape(x,[],1); % ensure x is a single column


hvdl = @haarhoffvdlinde;           fignum = 6231;    % choose: non-symbolic
% hvdl = @haarhoffvdlinde_symbolic;  fignum = 6232;    % choose: symbolic

% prepare figure
figure(fignum)
clf
y = hvdl(x,p);
plot(x,y,'.-');
drawnow

% Symbolische Ableitung auswerten
[yy, dydp, dydx] = hvdl(x,p); 
dyda0 = dydp(:,1);
dyda1 = dydp(:,2);
dyda2 = dydp(:,3);
dyda3 = dydp(:,4);

pause
clf 
disp('Showing analytical derivatives')
subplot(3,2,1); hold on; plot(x,     y, '-'); ylabel('y')
subplot(3,2,2); hold on; plot(x,  dydx, '-'); ylabel('dydx')
subplot(3,2,3); hold on; plot(x, dyda0, '-'); ylabel('dyda0')
subplot(3,2,4); hold on; plot(x, dyda1, '-'); ylabel('dyda1')
subplot(3,2,5); hold on; plot(x, dyda2, '-'); ylabel('dyda2')
subplot(3,2,6); hold on; plot(x, dyda3, '-'); ylabel('dyda3')
drawnow


% Finite Differenzen Approximation
FDh = 1.0d-8 * p;
FDdirs = eye(length(p));
y2a0 = hvdl(x,p+FDh(1)*FDdirs(:,1));
y2a1 = hvdl(x,p+FDh(2)*FDdirs(:,2));
y2a2 = hvdl(x,p+FDh(3)*FDdirs(:,3));
y2a3 = hvdl(x,p+FDh(4)*FDdirs(:,4));
dyda0_FD = (y2a0-y)/FDh(1);
dyda1_FD = (y2a1-y)/FDh(2);
dyda2_FD = (y2a2-y)/FDh(3);
dyda3_FD = (y2a3-y)/FDh(4);

FDx_h = 1e-5;
y2x  = hvdl(x+FDx_h, p);
dydx_FD = (y2x-y)/FDx_h;


pause
clf;
disp('Showing analytical derivatives and FD derivative approximations')
subplot(3,2,2); hold on; plot(x,  dydx, '-'); plot(x,  dydx_FD, '.'); ylabel('dydx')
subplot(3,2,3); hold on; plot(x, dyda0, '-'); plot(x, dyda0_FD, '.'); ylabel('dyda0')
subplot(3,2,4); hold on; plot(x, dyda1, '-'); plot(x, dyda1_FD, '.'); ylabel('dyda1')
subplot(3,2,5); hold on; plot(x, dyda2, '-'); plot(x, dyda2_FD, '.'); ylabel('dyda2')
subplot(3,2,6); hold on; plot(x, dyda3, '-'); plot(x, dyda3_FD, '.'); ylabel('dyda3')
drawnow


pause
clf
disp('Showing difference between analytical derivatives and FD approximations')
subplot(3,2,2); hold on; plot(x,  dydx-dydx_FD , '.-'); ylabel('dydx diff')
subplot(3,2,3); hold on; plot(x, dyda0-dyda0_FD, '.-'); ylabel('dyda0 diff')
subplot(3,2,4); hold on; plot(x, dyda1-dyda1_FD, '.-'); ylabel('dyda1 diff')
subplot(3,2,5); hold on; plot(x, dyda2-dyda2_FD, '.-'); ylabel('dyda2 diff')
subplot(3,2,6); hold on; plot(x, dyda3-dyda3_FD, '.-'); ylabel('dyda3 diff')
drawnow
