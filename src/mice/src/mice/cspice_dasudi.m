function cspice_dasudi( handle, first, last, data )
   switch nargin
      case 4
         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
         data   = zzmice_int(data);
      otherwise
         error ( 'Usage: cspice_dasudi( handle, first, last, data )' )
   end
   try
      mice('dasudi_c', handle, first, last, data);
   catch spiceerr
      rethrow(spiceerr)
   end