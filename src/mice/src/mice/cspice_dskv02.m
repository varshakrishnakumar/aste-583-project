function [vrtces] = cspice_dskv02( handle, dladsc, start, room )
   switch nargin
      case 4
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         start  = zzmice_int( start  );
         room   = zzmice_int( room, [1, int32(inf)/2] );
      otherwise
         error ( [ 'Usage: [vrtces] = ' ...
                   'cspice_dskv02( handle, dladsc, start, room ) ' ] )
   end
   try
      [vrtces] = mice( 'dskv02_c', ...
                       handle, dladsc, start, room );
   catch spiceerr
      rethrow(spiceerr)
   end