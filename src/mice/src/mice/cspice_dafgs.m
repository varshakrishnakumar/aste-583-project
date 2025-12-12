function [dc, ic] = cspice_dafgs( nd, ni)
   switch nargin
      case 2
         nd  = zzmice_int(nd);
         ni  = zzmice_int(ni);
      otherwise
         error ( 'Usage: [dc, ic] = cspice_dafgs( nd, ni)' )
   end
   try
      [dc, ic] = mice( 'dafgs_c', nd, ni );
   catch spiceerr
      rethrow(spiceerr)
   end