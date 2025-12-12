function [surf] = mice_srfxpt( method, target, et, abcorr, obsrvr, dref, dvec)
   switch nargin
      case 7
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         dref   = zzmice_str(dref);
         dvec   = zzmice_dp(dvec);
      otherwise
         error ( ['Usage: [_surf_ ] = '                                    ...
                  'mice_srfxpt( `method`, `target`,  '                     ...
                  '_et_, `abcorr`, `obsrvr`, '                             ...
                  '`dref`, dvec(6))']  )
   end
   try
      [surf] = mice( 'srfxpt_s', method, target, et,                       ...
                     abcorr,     obsrvr, dref,   dvec );
   catch spiceerr
      rethrow(spiceerr)
   end