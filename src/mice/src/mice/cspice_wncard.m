function [wncard] = cspice_wncard(window)
   switch nargin
      case 1
         window = zzmice_win(window);
      otherwise
         error ( 'Usage: [wncard] = cspice_wncard( window )' )
   end
   try
      [wncard] = mice('wncard_c', [zeros(6,1); window] );
   catch spiceerr
      rethrow(spiceerr)
   end