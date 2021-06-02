for i in {1..1001};
do
    sbatch <(cat submit1.sh | sed "s|ARGUMENT|$i|g" )
done



