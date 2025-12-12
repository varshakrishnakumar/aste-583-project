function cspice_clpool
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: cspice_clpool' )
   end
   try
      mice('clpool_c');
   catch spiceerr
      rethrow(spiceerr)
   end