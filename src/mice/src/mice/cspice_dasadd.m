function cspice_dasadd( handle, data )
   switch nargin
      case 2
         handle = zzmice_int(handle);
         data   = zzmice_dp(data);
      otherwise
         error ( 'Usage: cspice_dasadd( handle, data(n) )' )
   end
   try
      mice('dasadd_c', handle, data);
   catch spiceerr
      rethrow(spiceerr)
   end