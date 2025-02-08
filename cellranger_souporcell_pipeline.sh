cd '/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows'
mkdir TN1
cd TN1
mkdir rmd res data figs cds
cd ..

#login
ssh -i "/Users/skanaan/" -Y skanaan@rhino
## check for jobs
squeue -u skanaan

#change these variables depending on the run
export fq=/home/sfurlan/scratch/demux/220223_VH00738_19_AAAKMLLHV/AAAKMLLHV/outs/fastq_path/AAAKMLLHV
export id=TN1                                                                                              #CHANGE EVERY TIME
samps=(TN_BM TN_34)                                                                                 #CHANGE EVERY TIME
#ls -alh $fq #optionally check fastq folder


#change these depending on genome and protein usage
export transcriptome=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A
protein_run=true
export featurefile=/fh/fast/furlan_s/grp/refs/totalseq/TSCv2/features.csv
export PROTtag="_TSCL"

#change this depending on the sample tag and desired merge name for souporcell
export GEXtag="_5p__GEX"
export mergename="merge"

#prepare run
export wd="/fh/scratch/delete90/furlan_s/${id}"
mkdir $wd
cd $wd
#rm -R !("archive") #optionally clean folder

#setup environment
export PRIORPATH=$PATH
ml SAMtools
ml Python/3.8.2-GCCcore-9.3.0
export PATH=/home/sfurlan/software/cellranger-6.0.2:$PATH

#####CELLRANGER######
export ALLCELLRANGERRUNS=()
export ALLCELLRANGERJOBS=()
for samp in ${samps[@]}; do
if [ "$protein_run" = true ]; then
  feature_ref_line="ref,${featurefile}"
  samp_CSP="${samp}${PROTtag}"
  antibody_capture_line="${samp_CSP},${fq},any,${samp_CSP},Antibody Capture,"
else
  feature_ref_line=""
  antibody_capture_line=""
fi
samp_GEX="${samp}${GEXtag}"
export csv="${wd}/samp_${samp}.csv"
export samp
cat > $csv << EOL
[gene-expression]
ref,$transcriptome
expect-cells,10000
include-introns

[feature]
$feature_ref_line

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
$samp_GEX,$fq,any,$samp_GEX,Gene Expression,
$antibody_capture_line
EOL
RUNNUM=$(sbatch -n 1 -c 24 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='cellranger multi --id=$samp \
                   --csv=$csv --localcores=24 --localmem=480')
export ALLCELLRANGERRUNS+="$RUNNUM:"
done



#####SOUPORCELL######
#ALLCELLRANGERRUNS=$(sbatch --wrap='echo hello')
ALLCELLRANGERRUNS=${ALLCELLRANGERRUNS::-1}
LASTRUNTEMP="${ALLCELLRANGERRUNS//'Submitted batch job '}"
export ALLCELLRANGERJOBS="${LASTRUNTEMP//' on cluster gizmo'}"
echo $ALLCELLRANGERJOBS
export ALLBARCODERUNS=()
export ALLBARCODEJOBS=()
export OUT=${wd}/${mergename}
mkdir $OUT
cd $OUT
export K=1
labels=()
bamvar=()
bcvar=()
for sample in ${samps[@]}; do
bamvar+=($wd/$sample/outs/per_sample_outs/$sample/count/sample_alignments.bam)
RUNNUM=$(sbatch -n 1 -c 1 --dependency=afterok:$ALLCELLRANGERJOBS --wrap="cut -d ',' -f2 $wd/$sample/outs/per_sample_outs/$sample/count/sample_barcodes.csv | gzip > $wd/$sample/outs/per_sample_outs/$sample/count/barcodes.tsv.gz")
export ALLBARCODERUNS+="$RUNNUM:"
bcvar+=($wd/$sample/outs/per_sample_outs/$sample/count/barcodes.tsv.gz)
labels+=("${sample}_")
done
ALLBARCODERUNS=${ALLBARCODERUNS::-1}
LASTRUNTEMP="${ALLBARCODERUNS//'Submitted batch job '}"
export ALLBARCODEJOBS="${LASTRUNTEMP//' on cluster gizmo'}"
echo $ALLBARCODEJOBS
export bams=$(IFS=, ; echo "${bamvar[*]}")
export bcs=$(IFS=, ; echo "${bcvar[*]}")
export samples=$(IFS=, ; echo "${labels[*]}")
RUNNUM=$(sbatch -n 1 -c 1 --dependency=afterok:$ALLBARCODEJOBS -p campus-new -M gizmo --mem-per-cpu=40000MB --wrap='mergeBams \
    -i $bams \
    -l $samples \
    -b $bcs \
    -o $OUT')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
RUNNUM=$(sbatch -n 1 -c 24 -p campus-new -M gizmo --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='samtools sort -@ 24 out.bam -o out.sorted.bam')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
RUNNUM=$(sbatch -n 1 -c 12 -p campus-new -M gizmo --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='samtools index -@ 12 out.sorted.bam')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
export REF=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
export VCF=/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
export sif=/home/sfurlan/sifs/souporcell.sif
export SINGULARITY_BINDPATH="/fh/scratch,/fh/fast,/shared"
export OUTDIR=$OUT/souporcell_1
mkdir -p $OUTDIR
RUNNUM=$(sbatch -n 1 -c 35 -p campus-new -M gizmo  --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='/home/sfurlan/software/singularity/bin/singularity exec $sif souporcell_pipeline.py \
                  -i out.sorted.bam -b outbcs.tsv.gz -f $REF --common_variants $VCF --skip_remap True \
                  -t 35 -o $OUTDIR -k $K')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
ks=(1 2 3 4 5 6 7 8)
for k in ${ks[@]}; do
export K=$k
export OUTDIR=${OUT}/souporcell_${K}
mkdir -p $OUTDIR
RUNNUM=$(sbatch -n 1 -c 35 -p campus-new -M gizmo --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='/home/sfurlan/software/singularity/bin/singularity exec $sif /opt/souporcell/souporcell/target/release/souporcell -a souporcell_1/alt.mtx -r souporcell_1/ref.mtx -b outbcs.tsv.gz --min_alt 2 --min_ref 2 -k $K -t 35 > $OUTDIR/clusters_tmp.tsv 2> $OUTDIR/log.tsv')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN2="${LASTRUNTEMP//' on cluster gizmo'}"
sbatch -n 1 -c 1 -p campus-new -M gizmo  --dependency=afterok:$LASTRUN2 --mem-per-cpu=16000MB --wrap='/home/sfurlan/software/singularity/bin/singularity exec $sif troublet -a souporcell_1/alt.mtx -r souporcell_1/ref.mtx --clusters $OUTDIR/clusters_tmp.tsv > $OUTDIR/clusters.tsv'
done



