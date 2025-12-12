function [data] = cspice_dafgda( handle, baddr, eaddr)
   switch nargin
      case 3
         handle = zzmice_int(handle);
         baddr  = zzmice_int(baddr);
         eaddr  = zzmice_int(eaddr);
      otherwise
         error ( 'Usage: data = cspice_dafgda( handle, baddr, eaddr)' )
   end
   try
      [data] = mice( 'dafgda_c', handle, baddr, eaddr );
   catch spiceerr
      rethrow(spiceerr)
   end