function [trgepc, srfvec, phase, solar, emissn, visibl, lit] =             ...
                cspice_illum_plid_pl02( handle, dladsc, target,            ...
                                        et,     abcorr, obsrvr,            ...
                                        spoint, plid                )
   switch nargin
      case 8
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         target = zzmice_str( target );
         et     = zzmice_dp( et );
         abcorr = zzmice_str( abcorr );
         obsrvr = zzmice_str( obsrvr );
         spoint = zzmice_dp( spoint );
         plid   = zzmice_int( plid );
      otherwise
         error ( [ 'Usage: [trgepc, srfvec(3), phase, solar, emissn, '     ...
                   'visibl, lit] = cspice_illum_plid_pl02( handle, '       ...
                   ' dladsc(SPICE_DLA_DSCSIZ), `target`, '                 ...
                   'et, `abcorr`, `obsrvr`, spoint(3), plid )' ] )
   end
   try
      [ilumin, visibl, lit]  = mice( 'illum_plid_pl02_s',                  ...
                                     handle, dladsc, target,               ...
                                     et,     abcorr, obsrvr, spoint, plid);
      trgepc   = reshape( [ilumin(:).trgepc], 1, [] );
      srfvec   = reshape( [ilumin(:).srfvec], 3, [] );
      phase    = reshape( [ilumin(:).phase ], 1, [] );
      solar    = reshape( [ilumin(:).incdnc], 1, [] );
      emissn   = reshape( [ilumin(:).emissn], 1, [] );
      visibl   = zzmice_logical(visibl);
      lit      = zzmice_logical(lit);
   catch spiceerr
      rethrow(spiceerr)
   end