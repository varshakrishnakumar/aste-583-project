function [nv, np] = cspice_dskz02( handle, dladsc )
   switch nargin
      case 2
         handle = zzmice_int( handle );
      otherwise
         error ( [ 'Usage: [nv, np] = ' ...
                   'cspice_dskz02( handle, dladsc ) ' ] )
   end
   try
      [nv, np] = mice( 'dskz02_c', handle, dladsc  );
   catch spiceerr
      rethrow(spiceerr)
   end