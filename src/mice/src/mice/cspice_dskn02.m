function [normal] = cspice_dskn02( handle, dladsc, plid )
   switch nargin
      case 3
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         plid   = zzmice_int( plid   );
      otherwise
         error ( [ 'Usage: [normal(3)] = cspice_dskn02( ' ...
                                     'handle, dladsc, plid )' ]);
   end
   try
      [normal] = mice('dskn02_c', handle, dladsc, plid );
   catch spiceerr
      rethrow(spiceerr)
   end