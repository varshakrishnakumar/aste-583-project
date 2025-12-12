function [c] = cspice_wnintd( a, b )
   switch nargin
      case 2
         a    = zzmice_win(a);
         b    = zzmice_win(b);
      otherwise
         error ( 'Usage: [c] = cspice_wnintd( a, b )' )
   end
   try
      [c] = mice('wnintd_c', [zeros(6,1); a], [zeros(6,1); b] );
   catch spiceerr
      rethrow(spiceerr)
   end