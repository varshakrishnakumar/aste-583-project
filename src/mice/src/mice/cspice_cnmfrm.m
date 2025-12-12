function [frcode, frname, found] = cspice_cnmfrm(cname)
   switch nargin
      case 1
         cname = zzmice_str(cname);
      otherwise
         error ( ['Usage: [_frcode_, _`frname`_, _found_]' ...
                                             ' = cspice_cnmfrm(_`cname`_)'] )
   end
   try
      cnmfrm = mice( 'cnmfrm_s', cname ) ;
      frcode = reshape( [cnmfrm(:).code],  1, [] );
      frname = char( cnmfrm.name );
      found  = reshape( [cnmfrm(:).found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end