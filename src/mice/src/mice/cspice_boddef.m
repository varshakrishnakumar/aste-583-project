function cspice_boddef(name, code)
   switch nargin
      case 2
         name = zzmice_str(name);
         code = zzmice_int(code);
      otherwise
         error ( 'Usage: cspice_boddef(`name`, code)' )
   end
   try
      mice('boddef_c', name, code);
   catch spiceerr
      rethrow(spiceerr)
   end