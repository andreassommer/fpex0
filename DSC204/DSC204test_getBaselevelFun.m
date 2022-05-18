function blfun = DSC204test_getBaselevelFun(dscdata)
 % INPUT:   dscdata --> structure as returned by DSC204_readFile
 
 error('This does not work anymore due to revised blds format');
 
 % config
 fignum = 128;
 
 % load default detection settings
 blds = DSC204_getBaselineDetectionSettings();
 detectionArgs = {blds.reldevA, blds.reldevB, blds.reldevS2, blds.absdevA, blds.absdevB, blds.absdevS2};

 % quick accessors
 [t, T, uV, sf, mW] = DSC204_quickAccessors(dscdata, blds.Tmin, blds.Tmax);
 X = T;
 Y = uV;
 
 % detect linear ranges
 [~, peakPos] = max(Y);
 initlen_l = floor(blds.initfraction * peakPos);
 initlen_r = floor(blds.initfraction * (length(X)-peakPos));
 [idx_l, reg_l, initreg_l] = DSC204_detectLinearRange(X,Y,'left' ,initlen_l,detectionArgs{:});
 [idx_r, reg_r, initreg_r] = DSC204_detectLinearRange(X,Y,'right',initlen_r,detectionArgs{:});

 % get baselevel function
 bl_lin = DSC204_getBaselinePrimitive(X, Y, idx_l, reg_l.a, reg_l.b, idx_r, reg_r.a, reg_r.b, 'linear');
 bl_sig = DSC204_getBaselinePrimitive(X, Y, idx_l, reg_l.a, reg_l.b, idx_r, reg_r.a, reg_r.b, 'sigmoidal', 100);
 
 % plot
 figure(fignum); clf; 
 plot(X,Y,'b-',X,bl_sig(X),'r-',X,bl_lin(X),'c--')
 legend('Data', 'sigmoidal', 'linear');
 title(dscdata.fileSpec, 'Interpreter', 'none');
 
end