# @author: Vishal Koparde
# @modified: Skyler Kuhn (add output prefix) due to problem with '.' in output PATH
# @description: converts STAR aligned SE bam file to forward and reverse strand bigwigs

bam="$1"
outprefix="$2"
samtools view -H $1|grep "^@SQ"|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > ${bam}.genome

echo "samtools view -b -f 16 $bam | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome > ${bam}.fwd.bg && bedSort ${bam}.fwd.bg ${bam}.fwd.bg" > ${bam}.dobg
echo "samtools view -b -F 16 $bam | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome > ${bam}.rev.bg && bedSort ${bam}.rev.bg ${bam}.rev.bg" >> ${bam}.dobg
parallel --will-cite -j 2 < ${bam}.dobg

echo "bedGraphToBigWig ${bam}.fwd.bg ${bam}.genome ${outprefix}.fwd.bw" > ${bam}.dobw
echo "bedGraphToBigWig ${bam}.rev.bg ${bam}.genome ${outprefix}.rev.bw" >> ${bam}.dobw
parallel --will-cite -j 2 < ${bam}.dobw

rm -f ${bam}.fwd*.bg ${bam}.rev*.bg ${bam}.do* ${bam}.genome
