function [handle] = cspice_dlaopn( fname, ftype, ifname, ncomch )
   switch nargin
      case 4
         fname  = zzmice_str(fname);
         ftype  = zzmice_str(ftype);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int(ncomch);
      otherwise
         error ( [ 'Usage: [handle] = '                                     ...
                   'cspice_dlaopn( `fname`, `ftype`, `ifname`, ncomch )' ] )
   end
   try
      [handle] = mice('dlaopn_c', fname, ftype, ifname, ncomch);
   catch spiceerr
      rethrow(spiceerr)
   end