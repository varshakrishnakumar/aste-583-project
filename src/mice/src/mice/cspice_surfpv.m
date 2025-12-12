function [stx, found] = cspice_surfpv( stvrtx, stdir, a, b, c )
   switch nargin
      case 5
         stvrtx = zzmice_dp(stvrtx);
         stdir  = zzmice_dp(stdir);
         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);
      otherwise
         error ( [ 'Usage: [stx(6), found] = '                              ...
                   'cspice_surfpv( stvrtx(6), stdir(6), a, b, c )' ] )
   end
   try
      [stx, found] = mice('surfpv_c', stvrtx, stdir, a, b, c);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end