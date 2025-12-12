function [window_f] = cspice_wninsd( left, right, window_i )
   switch nargin
      case 2
         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
      case 3
         left     = zzmice_dp(left);
         right    = zzmice_dp(right);
         window_i = zzmice_win(window_i);
      otherwise
         error ( ['Usage: [window_f] = cspice_wninsd( left, right, ', ...
                                                     '[window_i] )'] )
   end
   if ( nargin == 2 )
      try
         [window_f] = mice( 'wninsd_c', left, right );
      catch
         rethrow(lasterror)
      end
   else
      try
         [window_f] = mice( 'wninsd_c', left, right );
      catch
         rethrow(lasterror)
      end
      window_f = cspice_wnunid( window_f, window_i );
   end