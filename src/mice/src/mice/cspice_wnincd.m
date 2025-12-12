function [wnincd] = cspice_wnincd( left, right, window )
   switch nargin
      case 3
         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
         window = zzmice_win(window);
      otherwise
         error( 'Usage: [wnincd] = cspice_wnincd( left, right, window )' )
      end
   try
      [wnincd] = mice( 'wnincd_c', left, right, [zeros(6,1); window] );
      [wnincd] = zzmice_logical(wnincd);
   catch spiceerr
      rethrow(spiceerr)
   end