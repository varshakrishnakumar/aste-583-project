function [srfpts, pltids] = cspice_llgrid_pl02( handle, dladsc, grid )
   switch nargin
      case 3
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         grid   = zzmice_dp( grid );
      otherwise
         error ( [ 'Usage: [srfpts(3), pltids] = '                         ...
            'cspice_llgrid_pl02( handle, dladsc(SPICE_DLA_DSCSIZ), grid )' ] )
   end
   try
      [srfpts, pltids] = mice( 'llgrid_pl02', handle, dladsc, grid );
   catch spiceerr
      rethrow(spiceerr)
   end