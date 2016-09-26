#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 26, 2015

Description:  Generates random synthetic reads from a reference genome.  Read names indicate their genomic origin.
Allows precise customization of things like insert size and synthetic mutation type, sizes, and rates.
Read names generated by this program are used by MakeRocCure (samtoroc.sh) and GradeSamFile (gradesam.sh).
They can also be used by BBMap (bbmap.sh) and BBMerge (bbmerge.sh) to automatically calculate
true and false positive rates, if the flag 'parsecustom' is used.

Usage:   randomreads.sh ref=<file> out=<file> length=<number> reads=<number>


Basic parameters:
out=null            Output file.  If blank, filename(s) will be autogenerated.
ref=null            Reference file.  Not needed if the reference is already indexed.
build=1             If multiple references are indexed in the same directory,
                    each needs a unique build ID.
midpad=300          Specifies space between scaffolds in packed index.
reads=0             Generate this many reads (or pairs).
minlength=100       Generate reads of up to this length.
maxlength=100       Generate reads of at least this length.
length=100          Generate reads of exactly this length.
overwrite=t         Set to false to disallow overwriting of existing files.
replacenoref=f      Set to true to replace Ns in the reference sequence 
                    with random letters.
simplenames=f       Set to true to generate read names that clearly indicate
                    genomic origin, without BBMap internal coordinates.
illuminanames=f     Set to true to have matching names for paired reads, 
                    rather than naming by location.
renamebyinsert=f    Insert the insert size into the name.
spaceslash=f        Set true to add a space before slash read pairnum.
prefix=null         Generated reads will start with this prefix, 
                    rather than naming by location.
seed=0              Use this to set the random number generator seed; 
                    use -1 for a random seed.

Pairing parameters:
paired=f            Set to true for paired reads.
interleaved=f       Set to true for interleaved output (rather than in two files).
mininsert=          Controls minimum insert length.  Default depends on read length.
maxinsert=          Controls maximum insert length.  Default depends on read length.
triangle=t          Make a triangular insert size distribution.
flat=f              Make a roughly flat insert size distribution..
superflat=f         Make a perfectly flat insert size distribution.
gaussian=f          Make a bell-shaped insert size distribution, with 
                    standard deviation of (maxinsert-mininsert)/6.
samestrand=f        Generate paired reads on the same strand.

Mutation parameters:
snprate=0           Add snps to reads with this probability (0-1).
insrate=0           Add insertions to reads with this probability (0-1).
delrate=0           Add deletions to reads with this probability (0-1).
subrate=0           Add contiguous substitutions to reads with this probability (0-1).
nrate=0             Add nocalls to reads with this probability (0-1).

Note: With a 'rate' of X, each read has an X chance of getting at least 
1 mutation, X^2 chance of 2+ mutations, X^3 chance of 3+ mutations, 
and so forth up to the maximum allowed number of mutations of that type.

maxsnps=3           Add at most this many snps per read.
maxinss=2           Add at most this many deletions per read.
maxdels=2           Add at most this many insertions per read.
maxsubs=2           Add at most this many contiguous substitutions per read.
maxns=0             Add at most this many blocks of Ns per read.

maxinslen=12        Max length of insertions.
maxdellen=400       Max length of deletions.
maxsublen=12        Max length of contiguous substitutions.
maxnlen=1           Min length of N blocks.

mininslen=1         Min length of insertions.
mindellen=1         Min length of deletions.
minsublen=2         Min length of contiguous substitutions.
minnlen=1           Min length of N blocks.

Illumina quality parameters:
maxq=36             Upper bound of quality values.
midq=32             Approximate average of quality values.
minq=28             Lower bound of quality values.
q=                  Sets maxq, midq, and minq to the same value.
adderrors=t         Add substitution errors based on quality values, 
                    after mutations.
qv=4                Vary the base quuality of reads by up to this much
                    to simulate tile effects.

PacBio quality parameters:
pacbio=f            Use a PacBio error model, rather than Illumina 
                    error model, and add PacBio errors after mutations.
pbmin=0.13          Minimum rate of PacBio errors for a read.
pbmax=0.17          Maximum rate of PacBio errors for a read.

Other Parameters:
overlap=1           Require reads to overlap scaffold end by at least this much.
banns=f             Do not generate reads over reference Ns.
randomscaffold=f    Choose random scaffolds without respect to length.
amp=1               Simulate highly-amplified MDA single-cell data by 
                    setting this to a higher number like 1000.
replacenoref=f      Replace intra- and inter-scaffold Ns with random bases.
#colorspace=f       Generate Solid colorspace reads.
pbadapter=          Add adapter sequence to some reads using this literal string.
fragadapter=        Add this sequence to paired reads with insert size 
                    shorter than read length.
fragadapter2=       Use this sequence for read 2.


Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the 
                    program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 
		    200 megs.  The max is typically 85% of physical memory.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
z2="-Xms1g"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

randomreads() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.7_64bit
		module load pigz
	fi
	local CMD="java -d64 $EA $z -cp $CP align2.RandomReads3 build=1 $@"
	echo $CMD >&2
	eval $CMD
}

randomreads "$@"
