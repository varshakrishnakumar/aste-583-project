function [pictur, ok, errmsg] = cspice_tpictr(sample)
   switch nargin
      case 1
         sample = zzmice_str(sample);
      otherwise
         error ( [ 'Usage: [`pictur`, ok, `errmsg`] = ', ...
                                        'cspice_tpictr( `sample` )' ] )
   end
   try
      [pictur, ok, errmsg] =  mice('tpictr_c', sample );
      ok = zzmice_logical(ok);
   catch spiceerr
      rethrow(spiceerr)
   end