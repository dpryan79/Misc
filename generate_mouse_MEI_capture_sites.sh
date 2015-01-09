#!/bin/bash
#fetch
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromOut.tar.gz

#Extract and clean up
tar xf chromOut.tar.gz
rm -f mm10.out

#Merge and filter repeats
for d in `ls -d */`
do
    sed 's/[()]//g' $d/*.out | grep -v "\-int" | awk '{if(NR>3) {
        if($9=="C") {
            start = $14
            left = $12
        } else {
            start = $12
            left = $14
        }
        len=$13+left
        cov=($13-start+1)/len
        if($11=="LINE/L1") {
            if(cov>=0.98 && $2<=2.0 && $4<=1.0 && $3<=1.0) print $0;
        } else if($11=="SINE/Alu") {
            if(cov>=0.98 && $2<=2.0 && $4<=1.0 && $3<=1.0) print $0;
        } else if($11=="LTR/ERV1" || $11=="LTR/ERVK" || $11=="LTR/ERVL-MaLR") {
            if(cov>=0.98 && $2<=1.0 && $4<=0.5 && $3<=0.5) print $0;
        }}}' >> mm10.out
    rm -rf $d
done

#Convert to RepeatMasker output to BED format
awk '{if($7-$6>800) {
        printf("%s\t%i\t%s\t%s:%s\n",$5,$6-1,$6+399,$11,$10)
        printf("%s\t%i\t%s\t%s:%s\n",$5,$7-400,$7,$11,$10)
    } else {
        printf("%s\t%i\t%s\t%s:%s\n",$5,$6-1,$7,$11,$10)
    }}' mm10.out > mm10.bed
