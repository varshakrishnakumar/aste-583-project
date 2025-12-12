function [nrec] = cspice_cknr02( handle, descr )
   switch nargin
      case 2
         handle = zzmice_int(handle);
         descr  = zzmice_dp(descr);
      otherwise
         error ( 'Usage: [nrec] = cspice_cknr02( handle, descr(5) )' )
   end
   try
      [nrec] = mice('cknr02_c', handle, descr);
   catch spiceerr
      rethrow(spiceerr)
   end