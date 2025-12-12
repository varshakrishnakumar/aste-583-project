function [result] = cspice_gftfov( inst,   target, tshape, tframe,         ...
                                   abcorr, obsrvr, step,   cnfine, nintvls )
   switch nargin
      case 9
         inst    = zzmice_str(inst);
         target  = zzmice_str(target);
         tshape  = zzmice_str(tshape);
         tframe  = zzmice_str(tframe);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         step    = zzmice_dp(step);
         cnfine  = zzmice_win(cnfine);
         nintvls    = zzmice_int(nintvls,    [1, int32(inf)/2] );
      otherwise
         error ( [ 'Usage: [result] = cspice_gftfov( `inst`, '             ...
                                   '`target`, `tshape`, `tframe` '         ...
                                   '`abcorr`, `obsrvr`, step, '            ...
                                   'cnfine, nintvls )' ] )
   end
   try
      [result] = mice('gftfov_c', inst,  target,  tshape,  tframe,         ...
                                  abcorr, obsrvr, step,                    ...
                                 [zeros(6,1); cnfine], nintvls);
   catch spiceerr
      rethrow(spiceerr)
   end