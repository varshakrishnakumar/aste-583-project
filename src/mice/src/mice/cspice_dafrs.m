function cspice_dafrs( dc, ic )
   switch nargin
      case 2
         dc  = zzmice_dp(dc);
         ic  = zzmice_int(ic);
      otherwise
         error ( 'Usage: cspice_dafrs( dc(nd), ic(ni) )' )
   end
   try
      mice('dafrs_c', dc, ic);
   catch spiceerr
      rethrow(spiceerr)
   end