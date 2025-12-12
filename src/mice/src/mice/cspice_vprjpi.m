function [vout, found] = cspice_vprjpi( vin, projpl, invpl )
   switch nargin
      case 3
         vin    = zzmice_dp( vin );
         projpl = zzmice_pln( projpl );
         invpl  = zzmice_pln( invpl );
      otherwise
         error ( ['Usage: [vout(3), found] = ' ...
                                  'cspice_vprjpi( vin(3), projpl, invpl )'] )
   end
   try
      [vout, found] = mice('vprjpi_c', vin, projpl, invpl );
   catch spiceerr
      rethrow(spiceerr)
   end