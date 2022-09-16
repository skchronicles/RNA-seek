# @author: Vishal Koparde
# @modified: Skyler Kuhn
# @description: converts STAR aligned PE bam file to forward and reverse strand bigwigs

bam="$1"
outprefix="$2"

samtools view -H $1|grep "^@SQ"|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > ${bam}.genome

echo "samtools view -b -f 128 -F 16 $bam | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome > ${bam}.fwd1.bg && bedSort ${bam}.fwd1.bg ${bam}.fwd1.bg" > ${bam}.dobg
echo "samtools view -b -f 80 $bam | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome > ${bam}.fwd2.bg && bedSort ${bam}.fwd2.bg ${bam}.fwd2.bg" >> ${bam}.dobg
echo "samtools view -b -f 144 $bam | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome > ${bam}.rev1.bg && bedSort ${bam}.rev1.bg ${bam}.rev1.bg" >> ${bam}.dobg
echo "samtools view -b -f 64 -F 16 $bam | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome > ${bam}.rev2.bg && bedSort ${bam}.rev2.bg ${bam}.rev2.bg" >> ${bam}.dobg
parallel --will-cite -j 4 < ${bam}.dobg

bedtools unionbedg -i ${bam}.fwd1.bg ${bam}.fwd2.bg|awk -F"\t" -v OFS="\t" '{sum=0;for (i=4;i<=NF;i+=1) {sum+=$i};print $1,$2,$3,sum}' - > ${bam}.fwd.bg
bedtools unionbedg -i ${bam}.rev1.bg ${bam}.rev2.bg|awk -F"\t" -v OFS="\t" '{sum=0;for (i=4;i<=NF;i+=1) {sum+=$i};print $1,$2,$3,sum}' - > ${bam}.rev.bg
bedSort ${bam}.fwd.bg ${bam}.fwd.bg
bedGraphToBigWig ${bam}.fwd.bg ${bam}.genome ${outprefix}.fwd.bw
bedSort ${bam}.rev.bg ${bam}.rev.bg
bedGraphToBigWig ${bam}.rev.bg ${bam}.genome ${outprefix}.rev.bw

rm -f ${bam}.fwd*.bg ${bam}.rev*.bg ${bam}.dob? ${bam}.genome
