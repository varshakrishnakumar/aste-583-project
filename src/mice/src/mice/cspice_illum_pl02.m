function [phase, solar, emissn] = cspice_illum_pl02( handle, dladsc,       ...
                                                     target, et,           ...
                                                     abcorr, obsrvr,       ...
                                                     spoint         )
   switch nargin
      case 7
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         target = zzmice_str( target );
         et     = zzmice_dp( et );
         abcorr = zzmice_str( abcorr );
         obsrvr = zzmice_str( obsrvr );
         spoint = zzmice_dp( spoint );
      otherwise
         error ( [ 'Usage: [phase, solar, emissn] = '                      ...
            'cspice_illum_pl02( handle, dladsc(SPICE_DLA_DSCSIZ), '        ...
                  '`target`, et, `abcorr`, `obsrvr`, spoint(3) )' ] )
   end
   try
      [phase, solar, emissn] = mice( 'illum_pl02',                         ...
                                     handle, dladsc, target,               ...
                                     et,     abcorr, obsrvr, spoint);
   catch spiceerr
      rethrow(spiceerr)
   end