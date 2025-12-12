function [ocltid] = cspice_occult ( targ1, shape1,   frame1, ...
                                    targ2, shape2,   frame2, ...
                                    abcorr,  obsrvr, et )
   switch nargin
      case 9
         targ1  = zzmice_str(targ1);
         shape1 = zzmice_str(shape1);
         frame1 = zzmice_str(frame1);
         targ2  = zzmice_str(targ2);
         shape2 = zzmice_str(shape2);
         frame2 = zzmice_str(frame2);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         et     = zzmice_dp(et);
      otherwise
         error ( ['Usage: [_output_state_] = ' ...
                  'cspice_occult( `targ1`, `shape1`, `frame1`,' ...
                  '`targ2`, `shape2`, `frame2`, `abcorr`, '...
                  '`obsrvr`, _et_)'] )
   end
   try
      [ocltid] = mice('occult_c', targ1, shape1,   frame1, ...
                                  targ2, shape2,   frame2, ...
                                  abcorr,  obsrvr, et );
   catch spiceerr
      rethrow(spiceerr)
   end