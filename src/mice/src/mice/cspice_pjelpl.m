function [elout] = cspice_pjelpl( elin, plane )
   switch nargin
      case 2
         elin  = zzmice_ell(elin);
         plane = zzmice_pln(plane);
      otherwise
         error ( 'Usage: [elout] = cspice_pjelpl( elin, plane )' )
   end
   try
      [elout] = mice( 'pjelpl_s', elin, plane );
   catch spiceerr
      rethrow(spiceerr)
   end