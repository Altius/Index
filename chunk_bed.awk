#!/usr/bin/env awk -f
BEGIN{n=1; nWritten=0; outf="chunk0001.bed"}{
  if(NR > 1){
    if($1 != prevChr || $2-prev3 > minChunkSep){
      if(nWritten >= minPerChunk){
        outf = sprintf("chunk%04d.bed", ++n);
        nWritten = 0;
      }
    }
    else{
      if($2-prev3 > minChunkSep2 && nWritten > maxPerChunk){
        outf = sprintf("chunk%04d.bed",++n);
        nWritten=0;
      }
    }
  }
  print $0 > outf;
  nWritten++; prevChr = $1; prev3 = $3;
}
