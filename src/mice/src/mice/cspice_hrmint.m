function [f, df] = cspice_hrmint( n, xvals, yvals, x )
   switch nargin
      case 4
         n     = zzmice_int(n);
         xvals = zzmice_dp(xvals);
         yvals = zzmice_dp(yvals);
         x     = zzmice_dp(x);
      otherwise
         error ( [ 'Usage: [f, df] = '                                     ...
                   'cspice_hrmint( n, xvals(n), yvals(n+1), x )' ] )
   end
   try
      [f, df] = mice('hrmint_c', n, xvals, yvals, x);
   catch spiceerr
      rethrow(spiceerr)
   end