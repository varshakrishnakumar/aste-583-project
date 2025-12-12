function cspice_pdpool( name, values )
   switch nargin
      case 2
         name   = zzmice_str(name);
         values = zzmice_dp(values);
      otherwise
         error ( 'Usage: cspice_pdpool( `name`, values(n) )' )
   end
   try
      mice('pdpool_c', name, values );
   catch spiceerr
      rethrow(spiceerr)
   end