function [et] = cspice_str2et(timstr)
   switch nargin
      case 1
         timstr = zzmice_str(timstr);
      otherwise
         error ( 'Usage: [_et_] = cspice_str2et(_`timstr`_)' )
   end
   try
      [et] = mice('str2et_c', timstr);
   catch spiceerr
      rethrow(spiceerr)
   end