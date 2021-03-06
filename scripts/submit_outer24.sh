#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=08:00:00
#PBS -N PAR_SQUARE_BENCH
#PBS -e error.square
#PBS -o all.outer.24

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=24
export OMP_STACKSIZE=64M

echo "DFS 24"
for i in `seq 0 8`;
do
    aprun -n 1 -d 24 -cc none ./build/matmul_bench_dfs -outer_prod_like $i
done

echo "HYBRID 24"
for i in `seq 1 8`;
do
    aprun -n 1 -d 24 -cc none ./build/matmul_bench_hybrid -outer_prod_like $i
done

#echo "BFS 24"
#for i in `seq 1 8`;
#do
#    aprun -n 1 -d 24 -cc none ./build/matmul_bench_bfs -outer_prod_like $i
#done
