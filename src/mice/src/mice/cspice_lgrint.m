function [lgrint] = cspice_lgrint( n, xvals, yvals, x )
   switch nargin
      case 4
         n     = zzmice_int(n, [1, int32(inf)]);
         xvals = zzmice_dp(xvals);
         yvals = zzmice_dp(yvals);
         x     = zzmice_dp(x);
      otherwise
         error ( [ 'Usage: [lgrint] = '                                    ...
                   'cspice_lgrint( n, xvals(n), yvals(n), x )' ] )
   end
   try
      [lgrint] = mice('lgrint_c', n, xvals, yvals, x);
   catch spiceerr
      rethrow(spiceerr)
   end