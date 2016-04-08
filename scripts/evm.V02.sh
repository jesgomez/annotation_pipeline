#!/bin/bash 

################################################################################################################################################################
##This script runs evm into a temp directory from the step of making partitions to the step of collecting the results. After that, it copies to the current
##directory two files with the outputs (.out and .gff3).
##
##Before calling the script make sure that  you have all this files in your current directory: genome.fa,proteins.gff3, transcripts.gff3,
##predictions.gff3 and a file (.txt) with the weights for evm.
##
##It needs to be called with one argument, the name of your weights file but without the extension (eg if your weitghts are in weights.txt you need to call ##the script like this):
##           evm.sh weights
##
## It needs to be run in a himem node
##
## Author: Jessica Gomez 
## Date: 13032014
#################################################################################################################################################################

DIR=$PWD;
cd $TMPDIR;
genome=$DIR/genome.fa;
proteins=$DIR/proteins.gff3;
transcripts=$DIR/transcripts.gff3;
predictions=$DIR/predictions.gff3;
version=$1;
weights=$version.txt;
weights_path=$DIR/$weights;

echo partitioning;
$EVM_PATH/EvmUtils/partition_EVM_inputs.pl --genome $genome --gene_predictions $predictions --protein_alignments $proteins --transcript_alignments $transcripts --segmentSize 2000000 --overlapSize 1000000 --partition_listing partitions_list.out;

echo writing commands;
$EVM_PATH/EvmUtils/write_EVM_commands.pl --genome $genome --weights $weights_path --gene_predictions $predictions --protein_alignments $proteins --transcript_alignments $transcripts --output_file_name evm_$version.out --partitions partitions_list.out > evm_$version.cmd;

/project/devel/aateam/bin/annotation_scripts/split_commands_file_evm.pl evm_$version.cmd 16;

echo running evm
for ((i=1; i<=16; i++))
  do
 $EVM_PATH/EvmUtils/execute_EVM_commands.pl evm_$version.$i.cmd &
 done;

wait

echo collecting outputs;
$EVM_PATH/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm_$version.out;

echo converting to gff3
$EVM_PATH/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm_$version.out --genome $genome;

lfs find . -name evm_$version.out | xargs cat > $DIR/evm_$version.out
lfs find . -name evm_$version.out.gff3 | xargs cat > $DIR/evm_$version.gff3
echo EVM done!
cd $DIR;

