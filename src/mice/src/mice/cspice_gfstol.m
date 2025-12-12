function cspice_gfstol( value )
  switch nargin
      case 1
         value = zzmice_dp(value);
      otherwise
         error ( 'Usage: cspice_gfstol(value)' )
   end
   try
      mice('gfstol_c', value);
   catch spiceerr
      rethrow(spiceerr)
   end