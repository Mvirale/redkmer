#!/bin/bash

REDKMERSTEP1=$(qsub 1_redkmer_QCandMitocleanup.bash)
echo $REDKMERSTEP1

REDKMERSTEP2=$(qsub -W depend=afterany:$REDKMERSTEP1 2_redkmer_pacbins.bash)
echo $REDKMERSTEP2

	REDKMERSTEP2R=$(qsub -W depend=afterany:$REDKMERSTEP2 2R_submission.bash)
	echo $REDKMERSTEP2R

REDKMERSTEP3=$(qsub -W depend=afterany:$REDKMERSTEP2R 3_redkmer_kmerGenerate.bash)
echo $REDKMERSTEP3

REDKMERSTEP4=$(qsub -W depend=afterany:$REDKMERSTEP3 4_redkmer_kmersbowtie.bash)
echo $REDKMERSTEP4

REDKMERSTEP5=$(qsub -W depend=afterany:$REDKMERSTEP4 5_redkmer_uniqueness.bash)
echo $REDKMERSTEP5

	REDKMERSTEP5R=$(qsub -W depend=afterany:$REDKMERSTEP5 5R_submission.bash)
	echo $REDKMERSTEP5R

REDKMERSTEP6=$(qsub -W depend=afterany:$REDKMERSTEP5R 6_redkmer_kmer2genome.bash)
echo $REDKMERSTEP6

	REDKMERSTEP6R=$(qsub -W depend=afterany:$REDKMERSTEP6 6R_submission.bash)
	echo $REDKMERSTEP6R


