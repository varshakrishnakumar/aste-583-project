function [c] = cspice_wndifd( a, b )
   switch nargin
      case 2
         a    = zzmice_win(a);
         b    = zzmice_win(b);
      otherwise
         error ( 'Usage: [c] = cspice_wndifd( a, b )' )
   end
   try
      [c] = mice('wndifd_c', [zeros(6,1); a], [zeros(6,1); b] );
   catch spiceerr
      rethrow(spiceerr)
   end