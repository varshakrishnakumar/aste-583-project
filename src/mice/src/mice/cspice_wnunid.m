function [c] = cspice_wnunid( a, b )
   switch nargin
      case 2
         a    = zzmice_win(a);
         b    = zzmice_win(b);
      otherwise
         error( 'Usage: [c] = cspice_wnunid( a, b )' )
   end
   try
      [c] = mice( 'wnunid_c', [zeros(6,1); a], [zeros(6,1); b] );
   catch spiceerr
      rethrow(spiceerr)
   end