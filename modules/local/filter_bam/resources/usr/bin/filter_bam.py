#! /usr/bin/env python3

import pysam
import sys

def main(bamfile, positions, doublestranded=False):
    #store all positions
    pos_check = []
    with open(positions, 'r') as pos_file:
        for line in pos_file:
            pos_check.append(int(line.replace('\n','')))

    #open the files
    infile = pysam.AlignmentFile(bamfile, 'rb')
    out3term = pysam.AlignmentFile('output.deaminated3.bam', 'wb', template=infile)

    #
    # main loop
    #
    revcomp_table = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

    # Iterate the bamfile
    for read in infile:
        # get read and ref sequences, mask ref-sequence if positions are matched
        seq = ''
        ref = ''
        for qpos, rpos, ref_base in read.get_aligned_pairs(with_seq=True):
            if qpos is None:
                seq += '-'   # deletion in read
            else:
                seq += read.query_sequence[qpos]
            if rpos is not None and (rpos + 1) in pos_check:
                # mask matched reference positions (1-based list)
                ref += 'N'
            else:
                # normal reference base or gap
                ref += str(ref_base).replace('None', '-')

        if read.is_reverse:
            ref = ref.translate(revcomp_table)[::-1]
            seq = read.get_forward_sequence()
        
        rlen = len(seq)
        # check if C>T (or G>A) substitutions in the first or last 3 postions
        if not doublestranded:
            mism = [
                n for n,x in enumerate(ref)
                if ((x, seq[n]) == ('c','T')) and (n<=2 or n>=rlen-3)
            ]
        else:
            mism = [
                n for n,x in enumerate(ref)
                if (
                    ((x, seq[n]) == ('c','T') and n<=2) or # check C>T in 5'
                    ((x, seq[n]) == ('g','A') and n>=rlen-3) # check G>A in 3'
                )
            ]

        # skip non-damaged reads
        if len(mism) == 0:
            continue

        #write damaged read to the file(s)
        out3term.write(read)

    infile.close()
    out3term.close()

if __name__ == "__main__":
    bamfile = sys.argv[1]
    positions = sys.argv[2]
    doublestranded = 'doublestranded' in sys.argv
    
    main(bamfile, positions=positions, doublestranded=doublestranded)