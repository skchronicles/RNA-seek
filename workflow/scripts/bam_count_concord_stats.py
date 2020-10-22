import sys
import pysam

infile = pysam.Samfile(sys.argv[1], "rb")
count = 0
keep = True
isproper_count = 0

for DNAread in infile.fetch():
    if DNAread.is_proper_pair:
        isproper_count += 1
    if DNAread.has_tag("NH"):
        tag_value=DNAread.get_tag("NH")
        if tag_value != 1:
            keep = False

    if keep:
        count+=1

    keep = True

infile.close()
print(isproper_count)
print(count)
