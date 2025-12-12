function cspice_dasadi( handle, data )
   switch nargin
      case 2
         handle = zzmice_int(handle);
         data   = zzmice_int(data);
      otherwise
         error ( 'Usage: cspice_dasadi( handle, data(n) )' )
   end
   try
      mice('dasadi_c', handle, data);
   catch spiceerr
      rethrow(spiceerr)
   end