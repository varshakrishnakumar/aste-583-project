function [n, found] = cspice_szpool( name )
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: [n, found] = cspice_szpool( `name` )' )
   end
   try
      [n, found] = mice('szpool_c', name);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end