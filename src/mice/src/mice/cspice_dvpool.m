function cspice_dvpool(name)
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: cspice_dvpool(_`name`_)' )
   end
   try
      mice('dvpool_c', name);
   catch spiceerr
      rethrow(spiceerr)
   end