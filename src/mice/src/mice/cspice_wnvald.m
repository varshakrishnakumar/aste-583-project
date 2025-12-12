function [window_f] = cspice_wnvald( window_i )
   switch nargin
      case 1
         window_i = zzmice_cell( window_i, 'double');
      otherwise
         error ( 'Usage: [window_f] = cspice_wnvald( window_i )' )
   end
   try
      [window_f] = mice('wnvald_c', window_i );
   catch spiceerr
      rethrow(spiceerr)
   end