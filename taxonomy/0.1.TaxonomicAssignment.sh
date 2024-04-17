#Here I assume you have blast (v2.12.0) & the nt database downloaded

#how to update the blast database with the perl script, these are already run so no need to do it again


#split fasta into chunks of 100 seqs for fast processing (do this for the 3 markers independently)
cd  ~/data/SeaChange_WP4/Metabarcoding/4.taxonomicAssignments
awk -v size=100 -v pre=EUK -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' EUK.cleaned.fasta
awk -v size=100 -v pre=MAM -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' MAM.cleaned.fasta
awk -v size=100 -v pre=RIZ -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' RIZ.cleaned.fasta



##EUK

for i in EUK.* ;do bname=$(basename $i); echo -e '#''!'"/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --time=00:60:00\n#SBATCH --mem-per-cpu=250G\ncd /datasets/globe_databases/ncbi/nt/20240215\nmodule load blast\nblastn -query ~/data/SeaChange_WP4/Metabarcoding/4.taxonomicAssignments/$bname -num_threads 1 -db nt -out ~/data/SeaChange_WP4/Metabarcoding/4.taxonomicAssignments/results/$bname.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200" > scripts/script.$bname.sh ;done

##RIZ
for i in RIZ.* ;do bname=$(basename $i); echo -e '#''!'"/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --time=00:60:00\n#SBATCH --mem-per-cpu=250G\ncd /datasets/globe_databases/ncbi/nt/20240215\nmodule load blast\nblastn -query ~/data/SeaChange_WP4/Metabarcoding/4.taxonomicAssignments/$bname -num_threads 1 -db nt -out ~/data/SeaChange_WP4/Metabarcoding/4.taxonomicAssignments/results/$bname.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200" > scripts/script.$bname.sh ;done

##MAM
for i in MAM.* ;do bname=$(basename $i); echo -e '#''!'"/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --time=00:60:00\n#SBATCH --mem-per-cpu=250G\ncd /datasets/globe_databases/ncbi/nt/20240215\nmodule load blast\nblastn -query ~/data/SeaChange_WP4/Metabarcoding/4.taxonomicAssignments/$bname -num_threads 1 -db nt -out ~/data/SeaChange_WP4/Metabarcoding/4.taxonomicAssignments/results/$bname.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200" > scripts/script.$bname.sh ;done

cd scripts
for i in script.* ;  do  sbatch $i;done


#combine output from all the split files to form one mega taxonomy file
cd ../results
cat EUK* > EUK.cleaned.rawtaxonomy.txt
cat RIZ* > RIZ.cleaned.rawtaxonomy.txt
cat MAM* > MAM.cleaned.rawtaxonomy.txt