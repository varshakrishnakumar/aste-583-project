function [found] = cspice_expool( name )
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: [found] = cspice_expool( `name` )' )
   end
   try
      [found] = mice('expool_c', name);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end