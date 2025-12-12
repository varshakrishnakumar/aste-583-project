function [window_f] = cspice_wnextd( side, window)
   switch nargin
      case 2
         side   = zzmice_str(side);
         window = zzmice_win(window);
      otherwise
         error ( 'Usage: [window_f] = cspice_wnextd( side, window )' )
   end
   try
      [window_f] = mice('wnextd_c', side, window );
   catch spiceerr
      rethrow(spiceerr)
   end