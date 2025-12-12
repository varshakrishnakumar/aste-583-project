function [f, df] = cspice_hrmesp( n, first, step, yvals, x )
   switch nargin
      case 5
         n     = zzmice_int(n);
         first = zzmice_dp(first);
         step  = zzmice_dp(step);
         yvals = zzmice_dp(yvals);
         x     = zzmice_dp(x);
      otherwise
         error ( [ 'Usage: [f, df] = '                                     ...
                   'cspice_hrmesp( n, first, step, yvals(n), x )' ] )
   end
   try
      [f, df] = mice('hrmesp_c', n, first, step, yvals, x);
   catch spiceerr
      rethrow(spiceerr)
   end