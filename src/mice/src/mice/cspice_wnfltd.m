function [window_f] = cspice_wnfltd( smlint, window_i )
   switch nargin
      case 2
         smlint   = zzmice_dp(smlint);
         window_i = zzmice_win(window_i);
      otherwise
         error ( 'Usage: [window_f] = cspice_wnfltd( smlint, window_i )' )
   end
   try
      [window_f] = mice( 'wnfltd_c', smlint, window_i );
   catch spiceerr
      rethrow(spiceerr)
   end