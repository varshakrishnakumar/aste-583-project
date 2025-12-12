function [nrec] = cspice_cknr03( handle, descr )
   switch nargin
      case 2
         handle = zzmice_int(handle);
         descr  = zzmice_dp(descr, true);
      otherwise
         error ( 'Usage: [nrec] = cspice_cknr03( handle, descr(5) )' )
   end
   try
      [nrec] = mice('cknr03_c', handle, descr);
   catch spiceerr
      rethrow(spiceerr)
   end