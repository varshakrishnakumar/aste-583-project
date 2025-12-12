function [result] = cspice_gfposc( target, frame, abcorr, obsrvr,  ...
                                   crdsys, coord, relate, refval,  ...
                                   adjust, step,  nintvls, cnfine )
   switch nargin
      case 12
         target  = zzmice_str(target);
         frame   = zzmice_str(frame);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         crdsys  = zzmice_str(crdsys);
         coord  = zzmice_str(coord);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);
      otherwise
         error ( [ 'Usage: [result] = cspice_gfposc( `target`, `frame`, '  ...
                                      '`abcorr`, `obsrvr`, `crdsys`, '     ...
                                      '`coord`, `relate`, refval, '        ...
                                      'adjust, step, nintvls, cnfine )' ] )
   end
   try
      [result] = mice('gfposc_c', target, frame, abcorr,  obsrvr,  ...
                                  crdsys, coord, relate,  refval,  ...
                                  adjust, step,  nintvls,          ...
                                  [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end