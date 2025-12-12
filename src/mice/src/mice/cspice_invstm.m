function [invmat] = cspice_invstm( mat )
   switch nargin
      case 1
         mat = zzmice_dp(mat);
      otherwise
         error ( 'Usage: [invmat(6,6)] = cspice_invstm( mat(6,6) )' )
   end
   try
      [invmat] = mice('invstm_c', mat);
   catch spiceerr
      rethrow(spiceerr)
   end