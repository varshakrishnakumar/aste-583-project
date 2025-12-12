function [sum] = cspice_dafps( dc, ic )
   switch nargin
      case 2
         dc = zzmice_dp(dc);
         ic = zzmice_int(ic);
      otherwise
         error ( 'Usage: [sum()] = cspice_dafps( dc(nd), ic(ni) )' )
   end
   try
      [sum] = mice('dafps_c', dc, ic);
   catch spiceerr
      rethrow(spiceerr)
   end