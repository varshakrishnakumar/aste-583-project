function cspice_furnsh(file)
   switch nargin
      case 1
         file = zzmice_str(file);
      otherwise
         error ( 'Usage: cspice_furnsh(_`file`_)' )
   end
   try
      mice('furnsh_c',file);
   catch spiceerr
      rethrow(spiceerr)
   end