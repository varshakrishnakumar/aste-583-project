function [nd, ni] = cspice_dafhsf( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( 'Usage: [nd, ni] = cspice_dafhsf( handle )' )
   end
   try
      [nd, ni] = mice('dafhsf_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end