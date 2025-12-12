function cspice_unload(file)
   switch nargin
      case 1
         file = zzmice_str( file );
      otherwise
         error ( 'Usage: cspice_unload(_`file`_)' )
   end
   try
      mice('unload_c', file)
   catch spiceerr
      rethrow(spiceerr)
   end