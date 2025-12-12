function [left, right] = cspice_wnfetd( window, n )
   switch nargin
      case 2
         window = zzmice_win(window);
         n      = zzmice_int(n);
      otherwise
         error( 'Usage: [left, right] = cspice_wnfetd( window, n )' )
      end
   try
      [left, right] = mice( 'wnfetd_c', [zeros(6,1); window], n );
   catch spiceerr
      rethrow(spiceerr)
   end