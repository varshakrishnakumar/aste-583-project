function [result] = cspice_gfpa( target,  illmn,  abcorr, obsrvr, ...
                                 relate,  refval, adjust, step,   ...
                                 nintvls, cnfine )
   switch nargin
      case 10
         target  = zzmice_str(target);
         illmn   = zzmice_str(illmn);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);
      otherwise
         error ( [ 'Usage: [result] = cspice_gfpa( `target`, `illmn`,  '    ...
                                                  '`abcorr, obsrvr`, '      ...
                                                  '`relate`, refval, '      ...
                                                  'adjust, step, nintvls, ' ...
                                                  'cnfine )' ] )
   end
   try
      [result] = mice('gfpa_c', target, illmn, abcorr, obsrvr, relate, ...
                                refval, adjust, step, nintvls,         ...
                                [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end