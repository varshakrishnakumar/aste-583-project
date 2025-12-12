function [spoint, trgepc, srfvec] = cspice_subpnt( method, target, et, ...
                                                   fixref, abcorr, obsrvr )
   switch nargin
      case 6
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
      otherwise
         error ( ['Usage: [_spoint_, _trgepc_, _srfvec_] = ' ...
                  'cspice_subpnt( `method`, `target`,'       ...
                  ' _et_, `fixref`, `abcorr`, `obsrvr`)']  )
   end
   try
      [subpnt] = mice('subpnt_s', method, target, et, fixref, abcorr, obsrvr);
      spoint   = reshape( [subpnt.spoint], 3, [] );
      trgepc   = reshape( [subpnt.trgepc], 1, [] );
      srfvec   = reshape( [subpnt.srfvec], 3, [] );
   catch spiceerr
      rethrow(spiceerr)
   end