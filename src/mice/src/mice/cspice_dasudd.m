function cspice_dasudd( handle, first, last, data )
   switch nargin
      case 4
         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
         data   = zzmice_dp(data);
      otherwise
         error ( 'Usage: cspice_dasudd( handle, first, last, data(N) )' )
   end
   try
      mice('dasudd_c', handle, first, last, data);
   catch spiceerr
      rethrow(spiceerr)
   end