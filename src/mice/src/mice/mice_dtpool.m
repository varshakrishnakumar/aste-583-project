function [value] = mice_dtpool(name)
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: [_value_] = mice_dtpool(_`name`_)' )
   end
   try
      [value] = mice('dtpool_s',name);
   catch spiceerr
      rethrow(spiceerr)
   end