function [handle] = cspice_cklpf( fname )
   switch nargin
      case 1
         fname = zzmice_str(fname);
      otherwise
         error ( 'Usage: [handle] = cspice_cklpf( `fname` )' )
   end
   try
      [handle] = mice('cklpf_c', fname);
   catch spiceerr
      rethrow(spiceerr)
   end