function DSC204test_detectLinearParts(dsc)
%
%
%  INPUT:  dsc --> DSC204-structure as returned by DSC204_readFile(s)
%
%
%

% config
fignum = 268;
Tmin = 70;
Tmax = inf;

absdevA  = inf; reldevA  = inf;
absdevB  = 0.2; reldevB  = inf;
absdevS2 = inf; reldevS2 = 0.5;
initfraction = 0.20;

% walk through data
for k = 1:length(dsc)

   % quick accessors
   [t, T, uV, sf, mW] = DSC204_quickAccessors(dsc(k), Tmin, Tmax);

   % select signal for the following stuff
   X = T;  
   Y = uV;
   
   % determine peak and ranges
   [~, peakPos] = max(Y);

   initlen_l = floor(initfraction * peakPos);
   initlen_r = floor(initfraction * (length(X)-peakPos));
   
   % detect linear ranges
   [stop_l, reg_l, initreg_l] = DSC204_detectLinearRange(X,Y,'left' ,initlen_l,reldevA,reldevB,reldevS2,absdevA,absdevB,absdevS2);
   [stop_r, reg_r, initreg_r] = DSC204_detectLinearRange(X,Y,'right',initlen_r,reldevA,reldevB,reldevS2,absdevA,absdevB,absdevS2);
   
   % plot data
   figure(fignum); clf; hold on;
   axh = gca();
   plot(X,Y);
   title(sprintf('Exp #%2d: %s  ---  Mass = %g mg  ---  Rate: %g %s', ...
                 k, dsc(k).fileSpec, dsc(k).mass, dsc(k).Tinfo.Tstep, dsc(k).Tinfo.Tstepunit), ...
         'Interpreter', 'none');
   
   % plot linear parts
   x = X(1:stop_l);
   y = reg_l.a + reg_l.b * x;
   plot(x,y,'r-',x(1),y(1),'r.',x(end),y(end),'r.')
   x = X(end:-1:stop_r); 
   y = reg_r.a + reg_r.b * x;
   plot(x,y,'r-',x(1),y(1),'r.',x(end),y(end),'r.')
   
   % statistic display
   showReg(axh, 'topleft', reg_l);
   showReg(axh, 'topright', reg_r);
   
   % pause
   fprintf('%s', dsc(k).fileSpec))
   fprintf('Press ESC to continue.\n')
   pause
   
end




   function h = showReg(axh, position, reg)
      % make the string
      fields = fieldnames(reg);
      str = '';
      for i=1:length(fields)
         valstr   = sprintf('%12.6f',reg.(fields{i}));
         fieldstr = sprintf('%7s', fields{i});
         str      = sprintf('%s\n%s: %s ', str, fieldstr, valstr);
      end
      str = sprintf('%s\n',str);

      % position it
      ax = axis(axh); 
      xmin = ax(1); xmax = ax(2); ymin = ax(3); ymax = ax(4);
      switch lower(position)
         case 'topleft'
            h = text(xmin, ymax, str, 'HorizontalAlignment','left' , 'VerticalAlignment','top','FontName','FixedWidth');
         case 'topright'
            h = text(xmax, ymax, str, 'HorizontalAlignment','right', 'VerticalAlignment','top','FontName','FixedWidth');
         case 'bottomleft'
            h = text(xmin, ymin, str, 'HorizontalAlignment','left' , 'VerticalAlignment','bottom','FontName','FixedWidth');
         case 'bottomright'
            h = text(xmax, ymin, str, 'HorizontalAlignment','right', 'VerticalAlignmentbottomtop','FontName','FixedWidth');
         otherwise
            error('Invalid position: %s', position);
      end
   end




end
 