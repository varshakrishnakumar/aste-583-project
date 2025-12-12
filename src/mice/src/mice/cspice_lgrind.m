function [p, dp] = cspice_lgrind( n, xvals, yvals, x )
   switch nargin
      case 4
         n     = zzmice_int(n, [1, int32(inf)]);
         xvals = zzmice_dp(xvals);
         yvals = zzmice_dp(yvals);
         x     = zzmice_dp(x);
      otherwise
         error ( [ 'Usage: [p, dp] = '                                     ...
                   'cspice_lgrind( n, xvals(n), yvals(n), x )' ] )
   end
   try
      [p, dp] = mice('lgrind_c', n, xvals, yvals, x);
   catch spiceerr
      rethrow(spiceerr)
   end