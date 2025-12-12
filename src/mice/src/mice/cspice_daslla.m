function [lastc, lastd, lasti] = cspice_daslla( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( 'Usage: [lastc, lastd, lasti] = cspice_daslla( handle )' )
   end
   try
      [lastc, lastd, lasti] = mice('daslla_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end