bedtools window -a ../Down.bed -b ../../../../AD.brain.TDP43-DRIP-ATAC-5hmC/CallPeakshTDP435hmC/MergePeakDEseq2/Down.bed >Down_Down.bed
bedtools window -a ../Down.bed -b ../../../../AD.brain.TDP43-DRIP-ATAC-5hmC/CallPeakshTDP435hmC/MergePeakDEseq2/Up.bed >Down_Up.bed
bedtools window -a ../Up.bed -b ../../../../AD.brain.TDP43-DRIP-ATAC-5hmC/CallPeakshTDP435hmC/MergePeakDEseq2/Down.bed >Up_Down.bed
bedtools window -a ../Up.bed -b ../../../../AD.brain.TDP43-DRIP-ATAC-5hmC/CallPeakshTDP435hmC/MergePeakDEseq2/Up.bed >Up_Up.bed

cut -f 1-3 Down_Down.bed >DRIP_Down_Down.bed
cut -f 1-3 Down_Up.bed >DRIP_Down_Up.bed
cut -f 1-3 Up_Down.bed >DRIP_Up_Down.bed
cut -f 1-3 Up_Up.bed >DRIP_Up_Up.bed

cut -f 4-6 Down_Down.bed >5hmC_Down_Down.bed
cut -f 4-6 Down_Up.bed >5hmC_Down_Up.bed
cut -f 4-6 Up_Down.bed >5hmC_Up_Down.bed
cut -f 4-6 Up_Up.bed >5hmC_Up_Up.bed
