#!/bin/bash

run_batch() {
local n=${1}_N${2}_T${3}_$4
local s=batch_$n.sh
cat << EOF > $s
#!/bin/bash
#SBATCH -N $2
#SBATCH --exclusive
#SBATCH --tasks-per-node=$3
#SBATCH -o $n.txt
#SBATCH -e stderr.txt
mpirun $1 $4
EOF
chmod +x $s
sbatch $s
}

for t in {,non-}blocking-ring; do
    for n in 10000000 20000000 40000000 80000000; do
        run_batch $t 1 4 $n
        run_batch $t 1 8 $n
        run_batch $t 2 8 $n
    done
done