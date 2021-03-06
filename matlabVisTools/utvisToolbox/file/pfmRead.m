function im = pfmRead( fname );
% im = pfmRead( fname )
%
% Load a pfm image into a MatLab matrix.  
% Compatability with Jepson pfm output format
%

[fid,msg] = fopen( fname, 'r' );
if (fid == -1)
  error(msg);
end

[pars type]= pnmReadHeader(fid);
if (pars==-1)
  fclose(fid);
  error([fname ': cannot parse pnm header']);
end

fileEndianness = 'UNKNOWN';
if (type == 'PL' | type == 'FP') 
  fileEndianness = 'LITTLE';
elseif (type == 'PB' | type == 'FU') 
  fileEndianness = 'BIG';
else
  fclose(fid);
  error('PFM file must be of type PL/PB');
end

xdim = pars(1);
ydim = pars(2);
%%% Maximum pixel value
maxval = pars(3);
%fprintf(1, 'original image size: cols %d rows %d\n', xdim, ydim)
sz = xdim * ydim;

count = 0; 
im = [];
if strcmp(endianness, fileEndianness)
  %% Reading file in same endianness format
  stat = fseek(fid, -sz*4, 'eof');
  if ~stat
    [im,count]  = fread(fid,sz,'float32');
  end
else
 %% Reading file in opposite endianness format

 %% Close file and reopen in the correct endianness format.
 fclose(fid);
 if strcmp(fileEndiannes, 'LITTLE')
   format = 'ieee-le';
 else % assume big endian
   format = 'ieee-be';
 end
 [fid,msg] = fopen( fname, 'r', format );
 if (fid < 0)
   error(msg);
 end

 %% Rewind to beginning of data
 stat = fseek(fid, -sz*4, 'eof');
 if ~stat
   [im,count]  = fread(fid,sz,'float32');
 end

end
fclose(fid);
	  
if (count == sz)
  im = reshape( im, xdim, ydim )';
else
  fprintf(1,'Warning: File ended %d words short: %s\n', sz-count, fname);
  im = reshape( [im ; zeros(sz-count,1)], xdim, ydim)';
end
