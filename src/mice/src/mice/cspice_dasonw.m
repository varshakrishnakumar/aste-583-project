function [handle] = cspice_dasonw( fname, ftype, ifname, ncomr )
   switch nargin
      case 4
         fname  = zzmice_str(fname);
         ftype  = zzmice_str(ftype);
         ifname = zzmice_str(ifname);
         ncomr  = zzmice_int(ncomr);
      otherwise
         error ( [ 'Usage: [handle] = '                                     ...
                   'cspice_dasonw( `fname`, `ftype`, `ifname`, ncomr )' ] )
   end
   try
      [handle] = mice('dasonw_c', fname, ftype, ifname, ncomr);
   catch spiceerr
      rethrow(spiceerr)
   end