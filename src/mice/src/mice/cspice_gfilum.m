function [result] = cspice_gfilum( method, angtyp,  target, illmn,  ...
                                   fixref, abcorr,  obsrvr, spoint, ...
                                   relate, refval,  adjust,         ...
                                   step,   nintvls, cnfine )
   switch nargin
      case 14
         method  = zzmice_str(method);
         angtyp  = zzmice_str(angtyp);
         target  = zzmice_str(target);
         illmn   = zzmice_str(illmn);
         fixref  = zzmice_str(fixref);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         spoint  = zzmice_dp(spoint);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);
      otherwise
         error ( [ 'Usage: [result] = cspice_gfilum( `method`, `angtyp`, '...
                                    '`target`, `illmn`, `fixref`, '       ...
                                    '`abcorr`, `obsrvr`, spoint[3], '     ...
                                    '`relate`, refval, adjust, '          ...
                                    'step, nintvls, cnfine, )' ] )
   end
   try
      [result] = mice('gfilum_c', method,  angtyp, target, illmn,  ...
                                  fixref,  abcorr, obsrvr, spoint, ...
                                  relate,  refval, adjust, step,   ...
                                  nintvls, [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end