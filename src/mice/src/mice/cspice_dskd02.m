function [values] = cspice_dskd02( handle, dladsc, item, start, room )
   switch nargin
      case 5
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         item   = zzmice_int( item   );
         start  = zzmice_int( start  );
         room   = zzmice_int( room, [1, int32(inf)/2] );
      otherwise
         error ( [ 'Usage: [values] = ' ...
                   'cspice_dskd02( handle, dladsc, item, start, room ) ' ] )
   end
   try
      [values] = mice( 'dskd02_c', ...
                       handle, dladsc, item, start, room );
   catch spiceerr
      rethrow(spiceerr)
   end