function [plane] = cspice_nvc2pl( normal, konst )
   switch nargin
      case 2
         normal = zzmice_dp(normal);
         konst  = zzmice_dp(konst);
      otherwise
         error ( ['Usage: [plane] = cspice_nvc2pl( normal(3), konst )'] )
   end
   try
      [plane] = mice('nvc2pl_c', normal, konst );
   catch spiceerr
      rethrow(spiceerr)
   end