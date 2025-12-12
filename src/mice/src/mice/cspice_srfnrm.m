function [normls] = cspice_srfnrm( method, target, et, fixref, srfpts )
   switch nargin
      case 5
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         srfpts = zzmice_dp(srfpts);
      otherwise
         error ( ['Usage: [normls] =  cspice_srfnrm( `method`, `target`, ' ...
                                                'et, `fixref`, srfpts )'] )
   end
   try
      [normls] = mice( 'srfnrm_c', method, target, et, fixref, srfpts );
   catch spiceerr
      rethrow(spiceerr)
   end