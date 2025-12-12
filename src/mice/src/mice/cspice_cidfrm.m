function [frcode, frname, found] = cspice_cidfrm( cent )
   switch nargin
      case 1
         cent = zzmice_int(cent);
      otherwise
         error ( ['Usage: [_frcode_, _`frname`_, _found_]' ...
                                             ' = cspice_cidfrm(_cent_)'] )
   end
   try
      cidfrm = mice( 'cidfrm_s', cent ) ;
      frcode = reshape( [cidfrm(:).code],  1, [] );
      frname = char( cidfrm.name );
      found  = reshape( [cidfrm(:).found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end