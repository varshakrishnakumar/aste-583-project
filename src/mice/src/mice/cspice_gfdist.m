function [result] = cspice_gfdist( target, abcorr, obsrvr, relate, refval, ...
                                   adjust, step, nintvls, cnfine )
   switch nargin
      case 9
         target  = zzmice_str(target);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);
      otherwise
         error ( [ 'Usage: [result] = cspice_gfdist( `target`, `abcorr`, '  ...
                                     '`obsrvr`, `relate`, refval, adjust, ' ...
                                     'step, nintvls, cnfine )' ] )
   end
   try
      [result] = mice('gfdist_c', target, abcorr, obsrvr, relate, ...
                                  refval, adjust, step, nintvls,  ...
                                  [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end