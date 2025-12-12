function [result] = cspice_gfsubc( target, fixref, method, abcorr, ...
                                   obsrvr, crdsys, coord, relate,  ...
                                   refval, adjust, step,  nintvls, ...
                                   cnfine )
   switch nargin
      case 13
         target  = zzmice_str(target);
         fixref  = zzmice_str(fixref);
         method  = zzmice_str(method);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         crdsys  = zzmice_str(crdsys);
         coord   = zzmice_str(coord);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);
      otherwise
         error ( [ 'Usage: [result] = cspice_gfsubc( `target`, `fixref`, ' ...
                                '`method`, `abcorr`, `obsrvr`, `crdsys`, ' ...
                                '`coord`, `relate`, refval, adjust, '      ...
                                'step, nintvls, cnfine )' ] )
   end
   try
      [result] = mice('gfsubc_c', target, fixref, method, abcorr, ...
                                  obsrvr, crdsys, coord, relate,  ...
                                  refval, adjust, step,  nintvls, ...
                                  [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end