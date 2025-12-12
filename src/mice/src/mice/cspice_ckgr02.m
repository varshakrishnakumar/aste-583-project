function [record] = cspice_ckgr02( handle, descr, recno )
   switch nargin
      case 3
         handle = zzmice_int(handle);
         descr  = zzmice_dp(descr);
         recno  = zzmice_int(recno);
      otherwise
         error ( [ 'Usage: [record(10)] = '                                ...
                   'cspice_ckgr02( handle, descr(5), recno )' ] )
   end
   try
      [record] = mice('ckgr02_c', handle, descr, recno);
   catch spiceerr
      rethrow(spiceerr)
   end