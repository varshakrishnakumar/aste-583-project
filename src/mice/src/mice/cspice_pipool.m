function cspice_pipool( name, ivals )
   switch nargin
      case 2
         name  = zzmice_str(name);
         ivals = zzmice_int(ivals);
      otherwise
         error ( 'Usage: cspice_pipool( `name`, ivals(n))' )
   end
   try
      mice('pipool_c', name, ivals );
   catch spiceerr
      rethrow(spiceerr)
   end