function [trgepc, obspos, lmbpts, pltids] =                                 ...
                                  cspice_limb_pl02( handle, dladsc, target, ...
                                                    et,     fixref, abcorr, ...
                                                    obsrvr, npts           )
   switch nargin
      case 8
         handle  = zzmice_int( handle );
         dladsc  = zzmice_int( dladsc );
         target  = zzmice_str( target );
         et      = zzmice_dp( et );
         fixref  = zzmice_str( fixref );
         abcorr  = zzmice_str( abcorr );
         obsrvr  = zzmice_str( obsrvr );
         npts    = zzmice_int( npts );
      otherwise
         error ( [ 'Usage: [trgepc, obspos(3), lmbpts(3,npts), '           ...
                   'pltids(npts)] = cspice_limb_pl02( handle, '            ...
                   'dladsc(SPICE_DLA_DSCSIZ), '                            ...
                   '`target`, et, `fixref`, `abcorr`, `obsrvr`, npts )' ] )
   end
   try
      [trgepc, obspos, lmbpts, pltids] = mice( 'limb_pl02',                ...
                                               handle, dladsc, target, et, ...
                                               fixref, abcorr, obsrvr, npts );
   catch spiceerr
      rethrow(spiceerr)
   end