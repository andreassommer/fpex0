function DSC204test_linreg()
   % Tests the linreg function.
   
   % setup
   fignum = 329;
   
   % coefficients and range
   a = 0.6;
   b = 1.3;
   c = 0.1;
   s2 = 1.26;
   Xmin  = 1;
   Xmax  = 30;
   Xstep = 0.001;
   fitfirst = 50;
   
   % linear and quadratic function
   y  = @(x) b*x + a;
   %y2 = @(x) a + b*x + c*x.^2;

   % first linear part, then quadratic
   y2breaks = [Xmin ; (Xmin+Xmax)/2 ; Xmax];
   y2coeffs = [0 b a+(y2breaks(1)*b) ; ...
               c b a+(y2breaks(2)*b)];
   y2pp = mkpp(y2breaks, y2coeffs);
   y2   = @(x) ppval(y2pp, x);
   
   
   % sample 
   X  = Xmin:Xstep:Xmax;
   Y  = y(X);
   Y2 = y2(X);
   ee = random('normal',0,sqrt(s2),size(X));
      
   
   % add error (same!)
   Ye  = Y  + ee;
   Y2e = Y2 + ee;
   
   % full regression
   regYfull  = DSC204_linreg(X,Ye);  % full regression
   regY2full = DSC204_linreg(X,Y2e);
      
   % plot Y
   figure(fignum); clf;
   subplot(2,2,1); 
   hY = gca;
   plot(hY,X,Ye,'r.', X,Y,'g-', X,regYfull.a+X*regYfull.b,'b-');
   showReg(hY, regYfull);
   axY = axis(hY);

   % plot Y2
   subplot(2,2,2)
   hY2 = gca;
   plot(hY2,X,Y2e,'r.', X,Y2,'g-', X,regY2full.a+X*regY2full.b,'b-');
   showReg(hY2, regY2full);
   axY2 = axis(hY2);
   
   % keep same axis limits in Y-plot
   axis(hY, axY2); axY = axY2;
   
   % partial regression with plot
   idx = 1:fitfirst;
   regYp  = DSC204_linreg(X(idx),Ye(idx));
   regY2p = DSC204_linreg(X(idx),Y2e(idx));
   for k = (idx(end)+1):length(X)
      regYp  = DSC204_linreg(X(k), Ye(k), regYp);
      regY2p = DSC204_linreg(X(k), Y2e(k), regY2p);
      
      subplot(2,2,3);
      plot(X(1:k),Ye(1:k),'r.',X(1:k),regYp.a+X(1:k)*regYp.b);
      axis(axY)
      showReg(gca, regYp);

      subplot(2,2,4);
      plot(X(1:k),Y2e(1:k),'r.',X(1:k),regY2p.a+X(1:k)*regY2p.b);
      axis(axY2)
      showReg(gca, regY2p);

      pause
   end
   
   
   % finito
   return
   
   
   
   function h = showReg(axh, reg)
      ax = axis(axh);
      fields = fieldnames(reg);
      str = '';
      for i=1:length(fields)
         valstr   = sprintf('%10.2f',reg.(fields{i}));
         fieldstr = sprintf('%10s', fields{i});
         str      = sprintf('%s\n%20s: %s ', str, fieldstr, valstr);
      end
      str = sprintf('%s\n',str);
      h  = text(ax(2), ax(3), str, 'HorizontalAlignment','right', 'VerticalAlignment','bottom','FontName','FixedWidth');
   end
   
end