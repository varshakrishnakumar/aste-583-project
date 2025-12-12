function [result] = cspice_wncomd( left, right, window)
   switch nargin
      case 3
         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
         window = zzmice_win(window);
      otherwise
         error ( 'Usage: [result] = cspice_wncomd( left, right, window)' )
   end
   try
      [result] = mice('wncomd_c', left, right, [zeros(6,1); window] );
   catch spiceerr
      rethrow(spiceerr)
   end