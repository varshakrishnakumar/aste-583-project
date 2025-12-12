function [record] = cspice_ckgr03( handle, descr, recno )
   switch nargin
      case 3
         handle = zzmice_int(handle);
         descr  = zzmice_dp(descr, true);
         recno  = zzmice_int(recno);
      otherwise
         error ( [ 'Usage: [record(8)] = '                                 ...
                   'cspice_ckgr03( handle, descr(5), recno )' ] )
   end
   try
      [record] = mice('ckgr03_c', handle, descr, recno);
   catch spiceerr
      rethrow(spiceerr)
   end