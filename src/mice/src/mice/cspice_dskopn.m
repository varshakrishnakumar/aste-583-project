function [handle] = cspice_dskopn( fname, ifname, ncomch )
   switch nargin
      case 3
         fname  = zzmice_str(fname);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int(ncomch);
      otherwise
         error ( ['Usage: [handle] = cspice_dskopn( `fname`, ',
                                          '`ifname`, ncomch )'] )
   end
   try
      [handle] = mice( 'dskopn_c', fname, ifname, ncomch );
   catch spiceerr
      rethrow(spiceerr)
   end