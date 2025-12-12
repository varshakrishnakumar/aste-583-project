function [dladsc, found] = cspice_dlabbs( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( [ 'Usage: [dladsc(SPICE_DLA_DSCSIZ), found] = '            ...
                   'cspice_dlabbs( handle )' ] )
   end
   try
      [dladsc, found] = mice('dlabbs_c', handle);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end