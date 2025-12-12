function cspice_kclear
   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: cspice_kclear' )
   end
   try
      mice('kclear_c');
   catch spiceerr
      rethrow(spiceerr)
   end