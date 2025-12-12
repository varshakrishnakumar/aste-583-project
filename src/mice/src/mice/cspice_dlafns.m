function [nxtdsc, found] = cspice_dlafns( handle, dladsc )
   switch nargin
      case 2
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
      otherwise
         error ( ['Usage: [nxtdsc(SPICE_DLA_DSCSIZ), found] = '            ...
                  'cspice_dlafns( handle, dladsc(SPICE_DLA_DSCSIZ) )'] )
   end
   try
      [nxtdsc, found] = mice('dlafns_c', handle, dladsc );
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end