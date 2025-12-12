function [prvdsc, found] = cspice_dlafps( handle, dladsc )
   switch nargin
      case 2
         handle = zzmice_int(handle);
         dladsc = zzmice_int(dladsc);
      otherwise
         error ( [ 'Usage: [prvdsc(SPICE_DLA_DSCSIZ), found] = '            ...
                   'cspice_dlafps( handle, dladsc(SPICE_DLA_DSCSIZ) )' ] )
   end
   try
      [prvdsc, found] = mice('dlafps_c', handle, dladsc);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end