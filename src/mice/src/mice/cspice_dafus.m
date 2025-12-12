function [dc, ic] = cspice_dafus( sum, nd, ni )
   switch nargin
      case 3
         sum = zzmice_dp(sum);
         nd  = zzmice_int(nd);
         ni  = zzmice_int(ni);
      otherwise
         error ( 'Usage: [dc(nd), ic(ni)] = cspice_dafus( sum(N), nd, ni )' )
   end
   try
      [dc, ic] = mice( 'dafus_c', sum, nd, ni );
   catch spiceerr
      rethrow(spiceerr)
   end