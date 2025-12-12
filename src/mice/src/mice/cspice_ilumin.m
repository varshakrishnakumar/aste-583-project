function [trgepc, srfvec, phase, incdnc, emissn] = cspice_ilumin( method, ...
                                                     target, et, fixref, ...
                                                     abcorr, obsrvr, spoint)
   switch nargin
      case 7
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         spoint = zzmice_dp(spoint);
      otherwise
         error ( [ 'Usage: [trgepc, srfvec(3), phase, incdnc, emissn] = ' ...
                               'cspice_ilumin( `method`, `target`, et, ' ...
                               '`fixref`, `abcorr`, `obsrvr`, spoint(3))' ] )
   end
   try
      [ilumin] = mice('ilumin_s', method, target, et, ...
                                  fixref, abcorr, obsrvr, spoint);
      trgepc   = reshape( [ilumin(:).trgepc], 1, [] );
      srfvec   = reshape( [ilumin(:).srfvec], 3, [] );
      phase    = reshape( [ilumin(:).phase ], 1, [] );
      incdnc   = reshape( [ilumin(:).incdnc], 1, [] );
      emissn   = reshape( [ilumin(:).emissn], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end