function [pnear, dist] = cspice_nplnpt( linpt, lindir, point )
   switch nargin
      case 3
         linpt  = zzmice_dp(linpt);
         lindir = zzmice_dp(lindir);
         point  = zzmice_dp(point);
      otherwise
         error ( ['Usage: [pnear(3), dist] = ' ...
                  'cspice_nplnpt( linpt(3), lindir(3), point(3) )'] )
   end
   try
      [nplnpt] = mice( 'nplnpt_s', linpt, lindir, point );
      pnear    = reshape( [nplnpt.pos], 3, [] );
      dist     = reshape( [nplnpt.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end