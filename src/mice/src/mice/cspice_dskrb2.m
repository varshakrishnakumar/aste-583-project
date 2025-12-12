function [mncor3, mxcor3] = cspice_dskrb2( vrtces, plates, corsys, corpar )
   switch nargin
      case 4
         vrtces = zzmice_dp(vrtces);
         plates = zzmice_int(plates);
         corsys = zzmice_int(corsys);
         corpar = zzmice_dp(corpar);
      otherwise
         error ( ['Usage: [mncor3, mxcor3] = '   ...
                  'cspice_dskrb2( vrtces(3,m), ' ...
                  'plates(3,n), corsys, corpar(SPICE_DSK_NSYPAR) ) '] )
   end
   try
      [mncor3, mxcor3] = mice( 'dskrb2_c', vrtces, plates, corsys, corpar );
   catch spiceerr
      rethrow(spiceerr)
   end