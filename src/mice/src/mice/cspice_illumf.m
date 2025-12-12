function [trgepc, srfvec, phase, incdnc, emissn, visibl, lit] = ...
             cspice_illumf( method, target, ilusrc, et, fixref, ...
                            abcorr, obsrvr, spoint)
   switch nargin
      case 8
         method = zzmice_str(method);
         target = zzmice_str(target);
         ilusrc = zzmice_str(ilusrc);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         spoint = zzmice_dp(spoint);
      otherwise
         error (                                                      ...
         [ 'Usage: '                                                  ...
         '[trgepc, srfvec(3), phase, incdnc, emissn, visibl, lit] = ' ...
         'cspice_ilumin( `method`, `target`, `ilusrc`, et, '          ...
         '`fixref`, `abcorr`, `obsrvr`, spoint(3))' ] )
   end
   try
      [ilumin, visibl, lit] = mice('illumf_s', method, target, ...
                                               ilusrc, et,     ...
                                               fixref, abcorr, ...
                                               obsrvr, spoint );
      trgepc   = reshape( [ilumin(:).trgepc], 1, [] );
      srfvec   = reshape( [ilumin(:).srfvec], 3, [] );
      phase    = reshape( [ilumin(:).phase ], 1, [] );
      incdnc   = reshape( [ilumin(:).incdnc], 1, [] );
      emissn   = reshape( [ilumin(:).emissn], 1, [] );
      visibl   = zzmice_logical(visibl);
      lit      = zzmice_logical(lit);
   catch spiceerr
      rethrow(spiceerr)
   end