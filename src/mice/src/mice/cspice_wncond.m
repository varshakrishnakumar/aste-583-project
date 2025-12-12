function [window_f] = cspice_wncond( left, right, window)
   switch nargin
      case 3
         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
         window = zzmice_win(window);
      otherwise
         error ( 'Usage: [window_f] = cspice_wncond( left, right, window)' )
   end
   try
      [window_f] = mice('wnexpd_c', -left, -right, window );
   catch spiceerr
      rethrow(spiceerr)
   end