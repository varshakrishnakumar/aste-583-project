function [trgsep] = cspice_trgsep( et, targ1, shape1, frame1, targ2, ...
                                   shape2, frame2, obsrvr, abcorr )
   switch nargin
      case 9
         et     = zzmice_dp(et);
         targ1  = zzmice_str(targ1);
         shape1 = zzmice_str(shape1);
         frame1 = zzmice_str(frame1);
         targ2  = zzmice_str(targ2);
         shape2 = zzmice_str(shape2);
         frame2 = zzmice_str(frame2);
         obsrvr = zzmice_str(obsrvr);
         abcorr = zzmice_str(abcorr);
      otherwise
         error ( [ 'Usage: [trgsep] = '                                     ...
                   'cspice_trgsep( et, `targ1`, `shape1`, `frame1`, '       ...
                   '`targ2`, `shape2`, `frame2`, `obsrvr`, `abcorr` )' ] )
   end
   try
      [trgsep] = mice('trgsep_c', et, targ1, shape1, frame1, targ2, shape2, ...
                      frame2, obsrvr, abcorr);
   catch spiceerr
      rethrow(spiceerr)
   end