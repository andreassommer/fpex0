% parameters
% a0 --> peak area
% a1 --> center of Gaussian component
% a2 --> standard deviation of the Gaussian component
% a3 --> peak distortion
p = [198 ; 120 ; 5 ; -1];
p = [228  ;  110  ;  3 ; -2.5];


% x range
x = 60:0.01:150;
x = reshape(x,[],1); % ensure x is a single column

% prepare figure
figure(6231)
clf
y = haarhoffvdlinde(x,p);
plot(x,y,'.-');

% Symbolische Ableitung auswerten
[yy, dydp] = haarhoffvdlinde(x,p); 
dyda0 = dydp(:,1);
dyda1 = dydp(:,2);
dyda2 = dydp(:,3);
dyda3 = dydp(:,4);

clf 
disp('Showing analytical derivatives')
subplot(3,2,1:2); plot(x,y); ylabel('y')
subplot(3,2,3); hold on; plot(x, dyda0, '-'); ylabel('dyda0')
subplot(3,2,4); hold on; plot(x, dyda1, '-'); ylabel('dyda1')
subplot(3,2,5); hold on; plot(x, dyda2, '-'); ylabel('dyda2')
subplot(3,2,6); hold on; plot(x, dyda3, '-'); ylabel('dyda3')
drawnow

return

% Finite Differenzen Approximation
FDh = 1.0d-4;
FDdirs = eye(length(p));
y2a0 = haarhoffvdlinde(x,p+FDh*FDdirs(:,1));
y2a1 = haarhoffvdlinde(x,p+FDh*FDdirs(:,2));
y2a2 = haarhoffvdlinde(x,p+FDh*FDdirs(:,3));
y2a3 = haarhoffvdlinde(x,p+FDh*FDdirs(:,4));

dyda0_FD = (y2a0-y)/FDh;
dyda1_FD = (y2a1-y)/FDh;
dyda2_FD = (y2a2-y)/FDh;
dyda3_FD = (y2a3-y)/FDh;

pause
clf;
disp('Showing analytical derivatives and FD derivative approximations')
subplot(3,2,1:2); plot(x,y); ylabel('y')
subplot(3,2,3); hold on; plot(x, dyda0, '-'); plot(x, dyda0_FD, '.'); ylabel('dyda0')
subplot(3,2,4); hold on; plot(x, dyda1, '-'); plot(x, dyda1_FD, '.'); ylabel('dyda1')
subplot(3,2,5); hold on; plot(x, dyda2, '-'); plot(x, dyda2_FD, '.'); ylabel('dyda2')
subplot(3,2,6); hold on; plot(x, dyda3, '-'); plot(x, dyda3_FD, '.'); ylabel('dyda3')
drawnow


pause
clf
disp('Showing difference between analytical derivatives and FD approximations')
subplot(2,2,1); hold on; plot(x, dyda0-dyda0_FD, '.-'); ylabel('dyda0')
subplot(2,2,2); hold on; plot(x, dyda1-dyda1_FD, '.-'); ylabel('dyda1')
subplot(2,2,3); hold on; plot(x, dyda2-dyda2_FD, '.-'); ylabel('dyda2')
subplot(2,2,4); hold on; plot(x, dyda3-dyda3_FD, '.-'); ylabel('dyda3')
drawnow
