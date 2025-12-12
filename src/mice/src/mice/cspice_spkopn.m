function [handle] = cspice_spkopn( fname, ifname, ncomch )
   switch nargin
      case 3
         fname  = zzmice_str(fname);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int(ncomch);
      otherwise
         error ( 'Usage: [handle] = cspice_spkopn(`fname`, `ifname`, ncomch)' )
   end
   try
      [handle] = mice( 'spkopn_c', fname, ifname, ncomch );
   catch spiceerr
      rethrow(spiceerr)
   end