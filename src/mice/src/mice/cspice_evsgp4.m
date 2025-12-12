function [state] = cspice_evsgp4( et, geophs, elems )
   switch nargin
      case 3
         et     = zzmice_dp(et);
         geophs = zzmice_dp(geophs);
         elems  = zzmice_dp(elems);
      otherwise
         error ( ['Usage: [_starg(6)_] = ' ...
                  'cspice_spkezr( et, _geophs(8)_, _elems(10_) )'] )
   end
   try
      [state] = mice( 'evsgp4_c', et, geophs, elems );
   catch spiceerr
      rethrow(spiceerr)
   end