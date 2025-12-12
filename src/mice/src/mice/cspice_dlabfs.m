function [dladsc, found] = cspice_dlabfs( handle )
   switch nargin
      case 1
         handle = zzmice_int( handle );
      otherwise
         error ( [ 'Usage: [dladsc(SPICE_DLA_DSCSIZ), found] = '            ...
                   'cspice_dlabfs( handle )' ] )
   end
   try
      [dladsc, found] = mice('dlabfs_c', handle );
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end