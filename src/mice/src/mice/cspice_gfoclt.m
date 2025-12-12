function [result] = cspice_gfoclt( occtyp, front,  fshape, fframe,         ...
                                   back,   bshape, bframe, abcorr,         ...
                                   obsrvr, step,   cnfine, nintvls )
   switch nargin
      case 12
         occtyp  = zzmice_str(occtyp);
         front   = zzmice_str(front);
         fshape  = zzmice_str(fshape);
         fframe  = zzmice_str(fframe);
         back    = zzmice_str(back);
         bshape  = zzmice_str(bshape);
         bframe  = zzmice_str(bframe);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         step    = zzmice_dp(step);
         cnfine  = zzmice_win(cnfine);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
      otherwise
         error ( [ 'Usage: [result] = cspice_gfoclt( `occtyp`, '           ...
                                           '`front`, `fshape`, `fframe`, ' ...
                                           '`back`, `bshape`, `bframe`, '  ...
                                           '`abcorr`, `obsrvr`, step, '    ...
                                           'cnfine, nintvls )' ] )
   end
   try
      [result] = mice('gfoclt_c',  occtyp, front,  fshape, fframe, ...
                                   back,   bshape, bframe,         ...
                                   abcorr, obsrvr,  step,          ...
                                   [zeros(6,1); cnfine], nintvls);
   catch spiceerr
      rethrow(spiceerr)
   end