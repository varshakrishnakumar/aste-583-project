function [result] = cspice_gfsep( targ1, shape1, frame1,          ...
                                  targ2, shape2, frame2,          ...
                                  abcorr, obsrvr, relate, refval, ...
                                  adjust, step, nintvls, cnfine )
   switch nargin
      case 14
         targ1   = zzmice_str(targ1);
         shape1  = zzmice_str(shape1);
         frame1  = zzmice_str(frame1);
         targ2   = zzmice_str(targ2);
         shape2  = zzmice_str(shape2);
         frame2  = zzmice_str(frame2);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);
      otherwise
         error ( [ 'Usage: [result] = cspice_gfsep( `targ1`, `shape1`, ' ...
                               '`frame1`, `targ2`, `shape2`, `frame2`, ' ...
                               '`abcorr`, `obsrvr`, `relate`, refval, '  ...
                               'adjust, step, nintvls, cnfine )' ] )
   end
   try
      [result] = mice('gfsep_c',  targ1, shape1, frame1,          ...
                                  targ2, shape2, frame2,          ...
                                  abcorr, obsrvr, relate, refval, ...
                                  adjust, step, nintvls,          ...
                                  [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end