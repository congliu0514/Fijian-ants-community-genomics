#!/bin/bash
#SBATCH --mail-user="liu.cong@oist.jp"
#SBATCH --job-name=examl_step1
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mem-per-cpu=20G
#SBATCH -p largemem
#SBATCH -c 2
#SBATCH -n 12
#SBATCH -t 60-20
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load examl/3.0.17
module load raxml/8.2.4

PHY=phylogeny_outgroup.phy
DIR="/home/l/liu-cong/work/fijian_rad_project/3-phylogeny/examl_new/bootstrap/tree/"
MODEL=GAMMA

NAME=$(echo `basename $PHY` | sed 's/.phy//')

mkdir $DIR/bsfiles
mkdir $DIR/starttrees

echo "Phylip file: $PHY"
echo "Output dir: $DIR"
echo "ExaML model: $MODEL"

echo "Generating BS Files ... "

mkdir $DIR/bsfiles/phy
mkdir $DIR/bsfiles/pstree
mkdir $DIR/bsfiles/binary

raxmlHPC -T 24 -N 100 -b $RANDOM -f j -m GTRCAT \
	-s $PHY -n REPS -w $DIR/bsfiles/phy > /dev/null

for i in {0..99}
do
	BSPHY=$DIR/bsfiles/phy/*".phy.BS"$i
	raxmlHPC -T 24 -y -p $RANDOM -m GTRCAT \
		-s $BSPHY -w $DIR/bsfiles/pstree -n "BS"$i > /dev/null
	cd $DIR/bsfiles/binary
	parse-examl -s $BSPHY -m DNA -n $i"_bin" > /dev/null
	cd -
done

echo "ExaML Bootstrap tree searches"
mkdir $DIR/bstree

for i in {0..99}
do
        echo "Computing tree $i"
        srun --mpi=pmi2 examl-AVX -m $MODEL -S -t $DIR/bsfiles/pstree/"RAxML_parsimonyTree.BS"$i \
                -s $DIR/bsfiles/binary/$i"_bin.binary" -n "bootstrapTree"$i \
                        -w $DIR/bsfiles/bstree > /dev/null
        rm bstree/ExaML_binaryCheckpoint*
done

echo ""
echo "ExaML Runs Complete 0-100"

cat $DIR/bsfiles/bstree/ExaML_result* > $DIR/ExaML_bootstrap.trees

raxmlSSE3 -T 24 -m GTRGAMMA -f b \
	-t $DIR/"ExaML_result."$NAME"_FullML" -z $DIR/ExaML_bootstrap.trees -n $NAME

echo "ExaML Runs all done"