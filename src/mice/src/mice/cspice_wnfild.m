function [window_f] = cspice_wnfild( smlgap, window_i)
   switch nargin
      case 2
         smlgap   = zzmice_dp(smlgap);
         window_i = zzmice_win(window_i);
      otherwise
         error ( 'Usage: [window_f] = cspice_wnfild( smlgap, window_i )' )
   end
   try
      [window_f] = mice('wnfild_c', smlgap, window_i );
   catch spiceerr
      rethrow(spiceerr)
   end