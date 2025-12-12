function [point, found] = cspice_surfpt ( positn, u, a, b, c )
   switch nargin
      case 5
         positn  = zzmice_dp(positn);
         u       = zzmice_dp(u);
         a       = zzmice_dp(a);
         b       = zzmice_dp(b);
         c       = zzmice_dp(c);
      otherwise
         error ( [ 'Usage: [point(3), found] = cspice_surfpt( positn(3), ' ...
                                                  'u(3), a, b, c )' ] )
   end
   try
      [surfpt] = mice('surfpt_s', positn, u, a, b, c );
      point  = reshape( [surfpt.spoint], 3, [] );
      found  = reshape( [surfpt.found ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end