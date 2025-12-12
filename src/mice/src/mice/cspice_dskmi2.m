function [spaixd, spaixi] = cspice_dskmi2( vrtces, plates, finscl, corscl, ...
                                           worksz, voxpsz, voxlsz, makvtl, ...
                                           spxisz )
   switch nargin
      case 9
         vrtces = zzmice_dp(vrtces);
         plates = zzmice_int(plates);
         finscl = zzmice_dp(finscl);
         corscl = zzmice_int(corscl);
         worksz = zzmice_int(worksz);
         voxpsz = zzmice_int(voxpsz);
         voxlsz = zzmice_int(voxlsz);
         makvtl = zzmice_int(makvtl);
         spxisz = zzmice_int(spxisz);
      otherwise
         error ( ['Usage: [ spaixd(SPICE_DSK02_IXDFIX), spaixi(spxisz)] = ' ...
                  'cspice_dskmi2( vrtces(3,m), plates(3,n), '               ...
                  'finscl, corscl, worksz, voxpsz, voxlsz, makvtl, spxisz)'] )
   end
   try
      [spaixd, spaixi] = mice( 'dskmi2_c', vrtces, plates, finscl, corscl, ...
                                           worksz, voxpsz, voxlsz, makvtl, ...
                                           spxisz);
   catch spiceerr
      rethrow(spiceerr)
   end