#!/usr/bin/env python
import os
import json
import argparse
import sys

#Author: Jessica Gomez, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
#Date:03242015



class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor"""
        #GENERAL PARAMETERS
        self.jsonFile = None                               #Name of the json configuration file to be created.
        self.run_geneid = True                             #By default run geneid step
        self.run_augustus = True                           #By default run augustus step
        self.run_genemark = True                           #By default run genemark step
        self.run_glimmer = True                            #By default run glimmer step
        self.run_geneid_introns = True                     #By default run geneid with introns step
        self.run_augustus_hints = True                   #By default run augustus with hints step
        self.run_genemark_ET = True                        #By default run genemark-ET step
        self.run_spaln = True                              #By default run spaln step
        self.run_pasa = True                               #By default run pasa step
        self.run_transdecoder = True                       #By default run transdecoder step
        self.run_evm = True                                #By default run EvidenceModeler step
        self.run_update = True                             #By default run annotation update step
        self.run_ncRNA_annotation = True                   #By default run the ncRNA annotation step          
        self.pipeline_HOME = "/project/devel/aateam/src/Annotation_pipeline.V03"   # Path to the pipeline home directory
 
        #INPUT PARAMETERS
        self.genome_masked = None                          #Path to the fasta genome masked.
        self.genome = None                                 #Path to the fasta genome.
        self.junctions = None                              #Path to the junctions gff file to run gene predictors with introns.
        self.incoding_junctions = None                     #Path to the junctions in coding regions gff file to run gene predictors with introns.
        self.species = None                                #Species name to run augustus with its trained parameters. For augustus and augustus with hints steps.
        self.geneid_parameters = None                      #Path to the geneid parameters file. For geneid, geneid with introns and framefixing (part of annotation update) steps.
        self.glimmer_directory = None                      #Path to the directory containing the trained parameters for running glimer.
        self.proteins = None                               #Path to the fasta with protein evidence.
        self.extrinsic_file_augustus_hints = None        #Extrinsic file to use when running augustus with hints. For more information read augustus documentation.
        self.ep_hints = None                               #Exonic hints to be used when running Augustus with hints. 
        self.pasadb = None                                 #Name of the pasa database, it must coincide with the name given in pasa-config.
        self.transcripts = None                            #Path to the fasta with transcript evidence.
        self.pasa_config = None                            #Path to the Pasa configuration file.
        self.RM_gff = None                                 #Path to the Repeat Masker gff output.
        self.cufflinks = None                              #Path to the cufflinks gtf file.
        self.update_config = None                          #Path to the Pasa configuration file for the annotation update step.
        self.project_name = None                           #Name of the project and version, to give the names to the final annotation output. 
        self.ncRNA_version = "None"                        #Version of the ncRNA annotation, to give the names to the final ncRNA annotation output. 

        #OUTPUT PARAMETERS
        self.annotation_step = "3"                         #Step of the annotation pipeline in the annotation process.    
        self.annotation_version = "01"                     #Version of the annotation process.      
        self.output_dir = "step0" + self.annotation_step + "_annotation_pipeline.V" + self.annotation_version  + "/"   #Directory to keep the outputs of the first annotation steps. 
        self.EVM_dir = "step0" + str(int(self.annotation_step) + 1) + "_EVM.V" + self.annotation_version  + "/"   #Directory to keep the files for the EVM step  
        self.dir_masked_chunks = self.output_dir + "/chunks_masked_reference"
        self.dir_genome_chunks = self.output_dir + "/chunks_genome_reference"
        self.dir_process_junctions = self.output_dir + "/junctions"         #Directory to keep all the files produced when getting incoding junctions.
        self.augustus_prediction = self.output_dir + "/gene_predictions/augustus/augustus_gene_prediction.gff3"  #Output file for the augustus predictions.
        self.augustus_preEVM = self.output_dir + "/gene_predictions/augustus/augustus_preEVM.gff3"  #Output file for the augustus predictions converted for EVM.
        self.geneid_prediction = self.output_dir + "/gene_predictions/geneid/geneid.gff3"  #Output file for the geneid predictions.
        self.geneid_preEVM = self.output_dir + "/gene_predictions/geneid/geneid_preEVM.gff3"  #Output file for the geneid predictions converted for EVM.
        self.genemark_prediction = self.output_dir + "/gene_predictions/genemark.gtf"  #Output file for the genemark predictions.
        self.genemark_preEVM = self.output_dir + "/gene_predictions/genemark_preEVM.gff3"  #Output file for the genemark predictions converted for EVM.
        self.glimmer_prediction = self.output_dir + "/gene_predictions/glimmer.gff3"  #Output file for the glimmer predictions.
        self.glimmer_preEVM = self.output_dir + "/gene_predictions/glimmer_preEVM.gff3"  #Output file for the glimmer predictions converted for EVM.
        self.spaln_cds =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_cds.gff3" #Output file for the spaln output in a cds gff3 format.
        self.spaln_gene =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_gene.gff3"    #Output file for the spaln output in a gene gff3 format.
        self.geneid_introns_prediction = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns.gff3"  #Output file for the geneid with introns predictions.
        self.geneid_introns_preEVM = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns_preEVM.gff3"  #Output file for the geneid with introns predictions converted for EVM.
        self.augustus_hints_prediction = self.output_dir + "/gene_predictions/augustus_with_hints/augustus_hints.gff3"  #Output file for the augustus with hints predictions.
        self.augustus_hints_preEVM = self.output_dir + "/gene_predictions/augustus_with_hints/augustus_hints_preEVM.gff3"  #Output file for the augustus with hints predictions converted for EVM.
        self.genemark_ET_prediction = self.output_dir + "/gene_predictions/genemark-ET.gtf"  #Output file for the genemark-ET predictions.
        self.genemark_ET_preEVM = self.output_dir + "/gene_predictions/genemark-ET_preEVM.gff3"  #Output file for the genemark-ET predictions converted for EVM.
        self.pasa_dir =  self.output_dir + "/protein_and_transcript_mappings/pasa/"     #Directory to keep all the pasa outputs.   
        self.update_dir = "step0" + str(int(self.annotation_step) + 2) + "_annotation_update.V" + self.annotation_version  + "/"   #Directory to keep the files for annotation update step.
        self.ncRNA_annotation_dir = "step0" + str(int(self.annotation_step) + 3) + "_ncRNA_annotation.V" + self.annotation_version  + "/"   #Directory to keep the files of the ncRNA annotation step.   
        self.out_cmsearch = self.ncRNA_annotation_dir + "/cmsearch.tbl"         #Output file to keep the cmsearch results
        self.out_tRNAscan = self.ncRNA_annotation_dir + "/tRNAscan-SE/tRNAscan.out"    # Output file to keep the tRNAscan-SE results. 

        #CHUNKS PARAMETERS
        self.masked_chunks = "50"                          #Number of chunks of the masked genome for parallelizing some gene predictors run.
        self.genome_chunks = "20"                          #Number of chunks of the genome for parallelizing tRNAscanSE.
        self.protein_chunks = "100"                        #Number of chunks to split the protein files for running blast and classify the lncRNAs. 

        #AUGUSTUS PARAMETERS
        self.aug_alternatives_from_sampling = "true"       #Report alternative transcripts generated through probabilistic sampling. For augustus and augustus with hints.                      
        self.aug_uniqueGeneId = "true"                     #If true, output gene identifyers like this: seqname.gN. For augustus and augustus with hints.         
        self.aug_gff3 = "ON"                               #Output in gff3 format. For augustus and augustus with hints.
        self.aug_sample = 60                               #For augustus and augustus with hints. 
        self.aug_noInFrameStop = "true"                    #Do not report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur. For augustus and augustus with hints.
        self.aug_maxtracks = 2                             #Maximum number of tracks allowed. For augustus and augustus with hints.
        self.aug_singlestrand = "false"                    #Predict genes independently on each strand, allow overlapping genes on opposite strands. For augustus and augustus with hints.
        self.aug_strand= "both"                            #For augustus and augustus with hints.
        self.aug_min_intron_len= 30                        #Minimum predicted intron length. For augustus and augustus with hints.                            
        self.augustus_weights = [2, 2, 1]                  #Weights given to augustus predictions when running EVM.
        self.additional_augustus_options = None            #Additional augustus options to run it, see augustus help for more information.

        #GENEID PARAMETERS
        self.geneid_weights = [2,1,2]                      #Weights given to geneid predictions when running EVM.
        self.geneid_options = "3U "                        #Desired geneid options to run it, see geneid documentation for more information.

       #GENEMARK PARAMETERS
        self.gmk_min_contig = 50000                        #Will ignore contigs shorter then min_contig in training
        self.gmk_max_contig = 5000000                      #will split input genomic sequence into contigs shorter than max_contig.
        self.gmk_max_gap = 5000                            #Will split sequence at gaps longer than max_gap. Letters 'n' and 'N' are       interpreted as standing within gaps 
        self.gmk_cores = 8                                 #Number of threads for running genemark. 
        self.additional_genemark_options = None            #Additional genemark options to run it, see genemark documentation for more information.
        self.genemark_weights = [1, 1, 1]                  #Weights given to genemark predictions when running EVM.

        #GLIMMER PARAMETERS
        self.glimmer_weights = [1, 1, 1]                   #Weights given to glimmer predictions when running EVM.

        #SPALN PARAMETERS
        self.spaln_ya = 1                                  #Stringency of splice site. 0->3: strong->weak.
        self.spaln_M = 4                                   #Number of outputs per query (1) (max=4 for genome vs cDNA|protein).
        self.spaln_O = 0                                   #0:Gff3_gene; 1:alignment; 2:Gff3_match; 3:Bed; 4:exon-inf; 5:intron-inf; 6:cDNA; 7:translated; 8: block-only; 12: binary.
        self.spaln_Q = 7                                   # 0:DP; 1-3: HSP-Search; 4-7; Block-Search.
        self.spaln_t = 8                                   #Number of threads.
        self.spaln_weights = [10, 8, 10]                   #Weights given to spaln mappings when running EVM.
        self.additional_spaln_options = None               #Additional spaln options to run it, see spaln help for more information.

        #GENEID INTRONS PARAMETERS
        self.geneid_introns_weights = [3,3,3]              #Weights given to geneid with intron predictions when running EVM.
        self.geneid_introns_options = "3nU"                #Desired geneid options to run geneid with introns, see geneid documentation for more information.

        #AUGUSTUS hints PARAMETERS
        self.augustus_hints_weights = [3,3,3]            #Weights given to augustus with intron predictions when running EVM.
        self.additional_augustus_hints_options = None    #Additional augustus options to run augustus with hints, see augustus help for more information.

        #GENEMARK-ET PARAMETERS
        self.gmk_et_score = 4                              #Minimum score of intron in initiation of the ET algorithm
        self.genemark_ET_weights = [3,3,3]                 #Weights given to augustus with intron predictions when running EVM.
        self.additional_genemark_ET_options = None         #Additional genemark-ET options to run it, see genemark documentation for more information.

        #PASA PARAMETERS
        self.pasa_CPU = 8                                  #Number of pasa_CPUs to run Pasa.
        self.pasa_step = 1                                 #Step from where to start running Pasa.
        self.pasa_weights = [8, 10, 8]                     #Weights given to pasa mappings when running EVM.
        self.pasa_home = "/project/devel/aateam/src/PASApipeline-2.0.2"
        self.create_database = False                       #By default do not create pasa database.

        #TRANSDECODER PARAMETERS
        self.transdecoder_weights = [3, 2, 3]              #Weights given to the pasa transdecoder output gff3 file.

        #EVM PARAMETERS
        self.evm_script = "/project/devel/aateam/src/Annotation_pipeline.V03/scripts/evm.V02.sh"

        #ncRNA ANNOTATION PARAMETERS
        self.cmsearch_CPUs = 16                            #Number of CPUs to run cmsearch
        self.Rfam = "/scratch/devel/talioto/de_novo_annotation/turbot/rna_annotation/CMs/Rfam.cm"    #CM file with the Rfam library.
###

        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.chunksParameters = {}
        self.augustusParameters = {}
        self.geneidParameters = {}
        self.genemarkParameters = {}
        self.glimmerParameters = {}
        self.spalnParameters = {}   
        self.geneidIntronsParameters = {}      
        self.augustusHintsParameters = {}
        self.genemarkETParameters = {}
        self.pasaParameters = {}
        self.transdecoderParameters = {}
        self.evmParameters = {}
        self.ncRNAannotationParameters = {}

####

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_chunks(parser)
        self.register_augustus(parser)
        self.register_geneid(parser)
        self.register_genemark(parser)
        self.register_glimmer(parser)
        self.register_spaln(parser)
        self.register_geneid_introns(parser)
        self.register_augustus_hints(parser)
        self.register_genemark_ET(parser)
        self.register_pasa(parser)
        self.register_transdecoder(parser)
        self.register_evm(parser)
        self.register_ncRNA_annotation(parser)

###

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('-jsonFile', dest="jsonFile", metavar="jsonFile", help='Configuration JSON to be generated. Default %s' % self.jsonFile)
        general_group.add_argument('--no-geneid', dest="run_geneid", action="store_false", help='If specified, do not run geneid step.')
        general_group.add_argument('--no-augustus', dest="run_augustus", action="store_false", help='If specified, do not run augustus step.')
        general_group.add_argument('--no-genemark', dest="run_genemark", action="store_false", help='If specified, do not run genemark step.')
        general_group.add_argument('--no-glimmer', dest="run_glimmer", action="store_false", help='If specified, do not run glimmer step.')
        general_group.add_argument('--no-geneid-introns', dest="run_geneid_introns", action="store_false", help='If specified, do not run geneid with introns step.')
        general_group.add_argument('--no-augustus-hints', dest="run_augustus_hints", action="store_false", help='If specified, do not run augustus with hints step.')
        general_group.add_argument('--no-genemark-ET', dest="run_genemark_ET", action="store_false", help='If specified, do not run genemark-ET step.')
        general_group.add_argument('--no-spaln', dest="run_spaln", action="store_false", help='If specified, do not run spaln step.')
        general_group.add_argument('--no-pasa', dest="run_pasa", action="store_false", help='If specified, do not run pasa step.')
        general_group.add_argument('--no-transdecoder', dest="run_transdecoder", action="store_false", help='If specified, do not run transdecoder step.')
        general_group.add_argument('--no-evm', dest="run_evm", action="store_false", help='If specified, do not run EVM step.')
        general_group.add_argument('--no-update', dest="run_update", action="store_false", help='If specified, do not run the annotation update step.')
        general_group.add_argument('--no-ncRNA-annotation', dest="run_ncRNA_annotation", action="store_false", help='If specified, do not run the ncRNA annotation step.')
        general_group.add_argument('--pipeline-HOME', dest="pipeline_HOME", metavar="pipeline_HOME", default='/project/devel/aateam/src/Annotation_pipeline.V03/', help='Path to the pipeline home directory. Default %s' %self.pipeline_HOME)

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--genome-masked', dest="genome_masked", metavar="genome_masked", help='Path to the fasta genome masked.')
        input_group.add_argument('--genome', dest="genome", metavar="genome", help='Path to the fasta genome.')
        input_group.add_argument('--junctions', dest="junctions", metavar="junctions", help='Path to the junctions gff file to run gene predictors with introns. Do not needed if incoding junctions existent file is given')
        input_group.add_argument('--incoding-junctions', dest="incoding_junctions", metavar="incoding_junctions", help='Path to the junctions in coding regions gff file to run gene predictors with introns. (Optional, it can be given if this file already exists or if you want to keep the processed incoding junctions in that concrete path.)')
        input_group.add_argument('--species', dest="species", metavar="species", help='Species name to run augustus with its trained parameters. For augustus and augustus with hints steps.')
        input_group.add_argument('--geneid-parameters', dest="geneid_parameters", metavar="geneid_parameters", help='Path to the geneid parameters file. For geneid, geneid with introns and framefixing (part of annotation update) steps.')
        input_group.add_argument('--glimmer-directory', dest="glimmer_directory", metavar="glimmer_directory", help='Path to the directory containing the trained parameters for running glimer.')
        input_group.add_argument('--proteins', dest="proteins", metavar="proteins", help='Path to the fasta with protein evidence.')
        input_group.add_argument('--extrinsic-file-augustus-hints', dest="extrinsic_file_augustus_hints", metavar="extrinsic_file_augustus_hints", help='Path to the Extrinsic file to use when running augustus with hints. For more information read augustus documentation.')
        input_group.add_argument('--ep-hints', dest="ep_hints", metavar="ep_hints", help='Path to the exonic hints to use when running augustus with hints. For more information read augustus documentation.')
        input_group.add_argument('--pasadb', dest="pasadb", metavar="pasadb", help='Name of the pasa database, it must coincide with the name in pasa_config.')
        input_group.add_argument('--transcripts', dest="transcripts", metavar="transcripts", help='Path to the fasta with transcript evidence.')
        input_group.add_argument('--pasa-config', dest="pasa_config", metavar="pasa_config", help='Path to the pasa configuration file.')
        input_group.add_argument('--RM-gff', dest="RM_gff", metavar="RM_gff", help='Path to the Repeat Masker gff output.')
        input_group.add_argument('--cufflinks', dest="cufflinks", metavar="cufflinks", help='Path to the cufflinks gtf file.')
        input_group.add_argument('--update-config', dest="update_config", metavar="update_config", help='Path to the Pasa configuration file.')
        input_group.add_argument('--project-name', dest="project_name", metavar="project_name", nargs=2, help='Name of the project and version of the annotation space separated, to give the names to the final annotation output.')
        input_group.add_argument('--ncRNA-version', dest="ncRNA_version", metavar="ncRNA_version", help='Version of the ncRNA annotation, to give the names to the final annotation output.')

    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--annotation-step', dest="annotation_step", default=3, type=int, help='Step of the annotation pipeline in the annotation process. Default %s' % self.annotation_step)
        output_group.add_argument('--annotation-version', dest="annotation_version", default='01', help='Version of the annotation process. Default %s' % self.annotation_version)
        output_group.add_argument('--output-dir', dest="output_dir", help='Directory to keep the outputs of the first annotation steps.  Default %s' % self.output_dir)
        output_group.add_argument('--EVM-dir', dest="EVM_dir", help='Directory to keep the files for the EVM step. Default %s' % self.EVM_dir)
        output_group.add_argument('--dir-masked-chunks', dest="dir_masked_chunks", help='Directory to keep the chunks of the masked genome. Default %s' % self.dir_masked_chunks)
        output_group.add_argument('--dir-genome-chunks', dest="dir_genome_chunks", help='Directory to keep the chunks of the genome. Default %s' % self.dir_genome_chunks)
        output_group.add_argument('--dir-process-junctions', dest="dir_process_junctions", help='Directory to keep all the files produced when getting incoding junctions. Default %s' % self.dir_process_junctions)
        output_group.add_argument('--augustus-prediction', dest="augustus_prediction",  help='Output file for the augustus predictions.  Default %s' % self.augustus_prediction)
        output_group.add_argument('--augustus-preEVM', dest="augustus_preEVM", help='Output file for the augustus predictions converted for EVM.  Default %s' % self.augustus_preEVM)
        output_group.add_argument('--geneid-prediction', dest="geneid_prediction", help='Output file for the geneid predictions.  Default %s' % self.geneid_prediction)
        output_group.add_argument('--geneid-preEVM', dest="geneid_preEVM", help='Output file for the geneid predictions converted for EVM.  Default %s' % self.geneid_preEVM)
        output_group.add_argument('--genemark-prediction', dest="genemark_prediction", help='Output file for the genemark predictions.  Default %s' % self.genemark_prediction)
        output_group.add_argument('--genemark-preEVM', dest="genemark_preEVM", help='Output file for the genemark predictions converted for EVM.  Default %s' % self.genemark_preEVM)
        output_group.add_argument('--glimmer-prediction', dest="glimmer_prediction", help='Output file for the glimmer predictions.  Default %s' % self.glimmer_prediction)
        output_group.add_argument('--glimmer-preEVM', dest="glimmer_preEVM", help='Output file for the glimmer predictions converted for EVM.  Default %s' % self.glimmer_preEVM)
        output_group.add_argument('--spaln-cds', dest="spaln_cds", help='Output file for the spaln output in a cds gff3 format.  Default %s' % self.spaln_cds)
        output_group.add_argument('--spaln-gene', dest="spaln_gene", help='Output file for the spaln output in a gene gff3 format.  Default %s' % self.spaln_gene)
        output_group.add_argument('--geneid-introns-prediction', dest="geneid_introns_prediction", help='Output file for the geneid with introns predictions.  Default %s' % self.geneid_introns_prediction)
        output_group.add_argument('--geneid-introns-preEVM', dest="geneid_introns_preEVM", help='Output file for the geneid with introns predictions converted for EVM.  Default %s' % self.geneid_introns_preEVM)
        output_group.add_argument('--augustus-hints-prediction', dest="augustus_hints_prediction", help='Output file for the augustus with hints predictions.  Default %s' % self.augustus_hints_prediction)
        output_group.add_argument('--augustus-hints-preEVM', dest="augustus_hints_preEVM", help='Output file for the augustus with hints predictions converted for EVM.  Default %s' % self.augustus_hints_preEVM)
        output_group.add_argument('--genemark-ET-prediction', dest="genemark_ET_prediction", help='Output file for the genemark-ET predictions.  Default %s' % self.genemark_ET_prediction)
        output_group.add_argument('--genemark-ET-preEVM', dest="genemark_ET_preEVM", help='Output file for the genemark-ET predictions converted for EVM.  Default %s' % self.genemark_ET_preEVM)
        output_group.add_argument('--pasa-dir', dest="pasa_dir", help='Directory to keep all the pasa outputs.  Default %s' % self.pasa_dir)
        output_group.add_argument('--update-dir', dest="update_dir", help='Directory to keep the files for annotation update step.   Default %s' % self.update_dir)
        output_group.add_argument('--ncRNA-annotation-dir', dest="ncRNA_annotation_dir", help='Directory to keep the files of the ncRNA annotation step. Default %s' % self.ncRNA_annotation_dir)
        output_group.add_argument('--out-cmsearch', dest="out_cmsearch", help='Output file to keep the cmsearch results. Default %s' % self.out_cmsearch)
        output_group.add_argument('--out-tRNAscan', dest="out_tRNAscan", help='Output file to keep the tRNAscan-SE results. Default %s' % self.out_tRNAscan)

    def register_chunks(self, parser):
        """Register all parameters for making chunks
        with the given argparse parser

        parser -- the argparse parser
        """
        chunks_group = parser.add_argument_group('Chunks')
        chunks_group.add_argument('--masked-chunks', dest="masked_chunks", type=int, default=50, help='Number of chunks of the masked genome for parallelizing some gene predictors run. Default %s' % self.masked_chunks)
        chunks_group.add_argument('--genome-chunks', dest="genome_chunks", type=int, default=20, help='Number of chunks of the genome for parallelizing cmsearch. Default %s' % self.genome_chunks)
        chunks_group.add_argument('--protein-chunks', dest="protein_chunks", type=int, default=100, help='Number of chunks to split the protein files for running blast and classify the lncRNAs.  Default %s' % self.protein_chunks)

    def register_augustus(self, parser):
        """Register all augustus parameters with given
        argparse parser

        parser -- the argparse parser
        """
        augustus_group = parser.add_argument_group('Augustus parameters')
        augustus_group.add_argument('--aug-alternatives-from-sampling', dest="aug_alternatives_from_sampling", default=self.aug_alternatives_from_sampling, choices=['true', 'false'], help='''Report alternative transcripts generated through probabilistic sampling. For augustus and augustus with hints. Default %s''' % str(self.aug_alternatives_from_sampling))
        augustus_group.add_argument('--aug-uniqueGeneId', dest="aug_uniqueGeneId", default=self.aug_uniqueGeneId, choices = ['true', 'false'], help='''If true, output gene identifyers like this: seqname.gN. For augustus and augustus with hints. Default %s''' % str(self.aug_uniqueGeneId))
        augustus_group.add_argument('--aug-gff3', dest="aug_gff3", default=self.aug_gff3, choices = ['ON', 'OFF', 'on', 'off'], help='''Output in gff3 format. For augustus and augustus with hints. Default %s''' % str(self.aug_gff3))
        augustus_group.add_argument('--aug-sample', dest="aug_sample", type=int, default=self.aug_sample, help='''For augustus and augustus with hints. Default %s''' % str(self.aug_sample))
        augustus_group.add_argument('--aug-noInFrameStop', dest="aug_noInFrameStop", default=self.aug_noInFrameStop, choices = ['true', 'false'], help='''Do not report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur. For augustus and augustus with hints. Default %s''' % str(self.aug_noInFrameStop))
        augustus_group.add_argument('--aug-maxtracks', dest="aug_maxtracks", type=int, default=self.aug_maxtracks, help='''Maximum number of tracks allowed. For augustus and augustus with hints. Default %s''' % str(self.aug_maxtracks))
        augustus_group.add_argument('--aug-singlestrand', dest="aug_singlestrand", default=self.aug_singlestrand, choices = ['true', 'false'], help='''Predict genes independently on each strand, allow overlapping genes on opposite strands. For augustus and augustus with hints. Default %s''' % str(self.aug_singlestrand))
        augustus_group.add_argument('--aug-strand', dest="aug_strand", default=self.aug_strand, choices=['both', 'forward', 'backward'], help='''For augustus and augustus with hints. Default %s''' % str(self.aug_strand))
        augustus_group.add_argument('--aug-min-intron-len', dest="aug_min_intron_len", type=int, default=self.aug_min_intron_len, help='''Minimum predicted intron length. For augustus and augustus with hints. Default %s''' % str(self.aug_min_intron_len))
        augustus_group.add_argument('--augustus-weights', dest="augustus_weights", nargs="+", type=int, default=self.augustus_weights, help='Weights given to augustus predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 2 2 1')
        augustus_group.add_argument('--additional-augustus-options', dest="additional_augustus_options", help='Additional augustus options to run it, see augustus help for more information about the possible options.')     

    def register_geneid(self, parser):
        """Register all geneid parameters with given
        argparse parser

        parser -- the argparse parser
        """
        geneid_group = parser.add_argument_group('Geneid parameters')
        geneid_group.add_argument('--geneid-weights', dest="geneid_weights", nargs="+", type=int, default=self.geneid_weights, help='Weights given to geneid predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 2 1 2 ')
        geneid_group.add_argument('--geneid-options', dest="geneid_options", default=self.geneid_options, help='Desired geneid options to run it, see geneid documentation for more information. Default %s''' % str(self.geneid_options))
  
    def register_genemark(self, parser):
        """Register all genemark parameters with given
        argparse parser

        parser -- the argparse parser
        """
        genemark_group = parser.add_argument_group('Genemark parameters')
        genemark_group.add_argument('--gmk-min-contig', dest="gmk_min_contig", type=int, default=self.gmk_min_contig, help='''Will ignore contigs shorter then min_contig in training. Default %s''' % str(self.gmk_min_contig))
        genemark_group.add_argument('--gmk-max-contig', dest="gmk_max_contig", type=int, default=self.gmk_max_contig, help='''Will split input genomic sequence into contigs shorter than max_contig. Default %s''' % str(self.gmk_max_contig))
        genemark_group.add_argument('--gmk-max-gap', dest="gmk_max_gap", type=int, default=self.gmk_max_gap, help='''Will split sequence at gaps longer than max_gap. Letters 'n' and 'N' are interpreted as standing within gaps. Default %s''' % str(self.gmk_max_gap))
        genemark_group.add_argument('--gmk-cores', dest="gmk_cores", type=int, default=self.gmk_cores, help='''Number of threads for running genemark. Default %s''' % str(self.gmk_cores))
        genemark_group.add_argument('--additional-genemark-options', dest="additional_genemark_options", help='Additional genemark options to run it, see genemark documentation for more information.')
        genemark_group.add_argument('--genemark-weights', dest="genemark_weights", nargs="+", type=int, default=self.genemark_weights, help='Weights given to genemark predictions when running EVM. Specify the weight for each EVM run separated by a space. Example  1 1 1 ') 

    def register_glimmer(self, parser):
        """Register all glimmer parameters with given
        argparse parser

        parser -- the argparse parser
        """
        glimmer_group = parser.add_argument_group('Glimmer parameters')
        glimmer_group.add_argument('--glimmer-weights', dest="glimmer_weights", nargs="+", type=int, default=self.glimmer_weights, help='Weights given to glimmer predictions when running EVM. Specify the weight for each EVM run separated by a space. Example  1 1 1 ') 

    def register_spaln(self, parser):
        """Register all spaln parameters with given
        argparse parser

        parser -- the argparse parser
        """
        spaln_group = parser.add_argument_group('Spaln parameters')
        spaln_group.add_argument('--spaln-ya', dest="spaln_ya", type=int, default = self.spaln_ya, choices = [0, 1, 2, 3], help='''Stringency of splice site. 0->3: strong->weak. Default %s''' % str(self.spaln_ya))
        spaln_group.add_argument('--spaln-M', dest="spaln_M", type=int, default = self.spaln_M, choices = [0, 1, 2, 3, 4], help='''Number of outputs per query (max=4 for genome vs cDNA|protein). Default %s''' % str(self.spaln_M))
        spaln_group.add_argument('--spaln-O', dest="spaln_O", type=int, default = self.spaln_O, choices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], help='''0:Gff3_gene; 1:alignment; 2:Gff3_match; 3:Bed; 4:exon-inf; 5:intron-inf; 6:cDNA; 7:translated; 8: block-only; 12: binary. Default %s''' % str(self.spaln_O))
        spaln_group.add_argument('--spaln-Q', dest="spaln_Q", type=int, default = self.spaln_Q, choices = [0, 1, 2, 3, 4, 5, 6, 7 ], help='''0:DP; 1-3: HSP-Search; 4-7; Block-Search. Default %s''' % str(self.spaln_Q))
        spaln_group.add_argument('--spaln-t', dest="spaln_t", type=int, default = self.spaln_t, help='''Number of threads. Default %s''' % str(self.spaln_t))
        spaln_group.add_argument('--spaln-weights', dest="spaln_weights", nargs="+", type=int, default=self.spaln_weights, help='Weights given to spaln mappings when running EVM. Specify the weight for each EVM run separated by a space. Example  10 8 10 ') 
        spaln_group.add_argument('--additional-spaln-options', dest="additional_spaln_options", help='Additional spaln options to run it, see spaln help for more information about the possible options.')

    def register_geneid_introns(self, parser):
        """Register all geneid with introns parameters with given
        argparse parser

        parser -- the argparse parser
        """
        geneid_introns_group = parser.add_argument_group('Geneid Introns parameters')
        geneid_introns_group.add_argument('--geneid-introns-weights', dest="geneid_introns_weights", nargs="+", type=int, default=self.geneid_introns_weights, help='Weights given to geneid with intron predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 3 3 3 ')
        geneid_introns_group.add_argument('--geneid-introns-options', dest="geneid_introns_options", default=self.geneid_introns_options, help='Desired geneid with intron options to run it, see geneid documentation for more information. Default %s''' % str(self.geneid_options))
  
    def register_augustus_hints(self, parser):
        """Register all augustus with hints parameters with given
        argparse parser

        parser -- the argparse parser
        """
        augustus_hints_group = parser.add_argument_group('Augustus hints parameters')
        augustus_hints_group.add_argument('--augustus-hints-weights', dest="augustus_hints_weights", nargs="+", type=int, default=self.augustus_hints_weights, help='Weights given to augustus with intron predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 3 3 3 ')
        augustus_hints_group.add_argument('--additional-augustus-hints-options', dest="additional_augustus_hints_options", default=self.additional_augustus_hints_options, help='Desired augustus with intron options to run it, see augustus documentation for more information.''')

    def register_genemark_ET(self, parser):
        """Register all genemark-ET parameters with given
        argparse parser

        parser -- the argparse parser
        """
        genemark_ET_group = parser.add_argument_group('Genemark-ET parameters')
        genemark_ET_group.add_argument('--gmk-et-score', dest="gmk_et_score", type=int, default=self.gmk_et_score, help='Minimum score of intron in initiation of the ET algorithm. Default %s''' %str(self.gmk_et_score)) 
        genemark_ET_group.add_argument('--additional-genemark-ET-options', dest="additional_genemark_ET_options", default=self.additional_genemark_ET_options, help='Additional genemark-ET options to run it, see genemark documentation for more information.')
        genemark_ET_group.add_argument('--genemark-ET-weights', dest="genemark_ET_weights", nargs="+", type=int, default=self.genemark_ET_weights, help='Weights given to genemark-ET predictions when running EVM. Specify the weight for each EVM run separated by a space. Example  3 3 3 ') 

    def register_pasa(self, parser):
        """Register all pasa parameters with given
        argparse parser

        parser -- the argparse parser
        """
        pasa_group = parser.add_argument_group('Pasa parameters')
        pasa_group.add_argument('--pasa-CPU', dest="pasa_CPU", type=int, default=self.pasa_CPU, help='''Number of pasa_CPUs to run Pasa. Default %s''' % str(self.pasa_CPU))
        pasa_group.add_argument('--pasa-step', dest="pasa_step", type=int, default=self.pasa_step, help='''Step from where to start running Pasa. Default %s''' % str(self.pasa_step))
        pasa_group.add_argument('--pasa-weights', dest="pasa_weights", nargs="+", type=int, default=self.pasa_weights, help='Weights given to pasa mappings when running EVM. Specify the weight for each EVM run separated by a space. Example  8 10 8 ')
        pasa_group.add_argument('--pasa-home', dest="pasa_home", default=self.pasa_home, help='''Path to the PASAHOME directory. Default %s''' %str(self.pasa_home ))
        pasa_group.add_argument('--create-database', dest="create_database", action="store_true", help='If specified, create pasa database.')

    def register_transdecoder(self, parser):
        """Register all transdecoder parameters with given
        argparse parser

        parser -- the argparse parser
        """
        transdecoder_group = parser.add_argument_group('Transdecoder parameters')
        transdecoder_group.add_argument('--transdecoder-weights', dest="transdecoder_weights", nargs="+", type=int, default=self.transdecoder_weights, help='Weights given to pasa transdecodergff3 output file  when running EVM. Specify the weight for each EVM run separated by a space. Example 3 2 3')

    def register_evm(self, parser):
        """Register all evm parameters with given
        argparse parser

        parser -- the argparse parser
        """
        evm_group = parser.add_argument_group('Evm parameters')
        evm_group.add_argument('--evm-script', dest="evm_script", default = self.evm_script, help='Script to run EVM. Default %s''' % str(self.evm_script))

    def register_ncRNA_annotation(self, parser):
        """Register all ncRNA annotation parameters with given
        argparse parser

        parser -- the argparse parser
        """
        ncRNA_group = parser.add_argument_group('ncRNA Annotation parameters')
        ncRNA_group.add_argument('--Rfam', dest="Rfam", default = self.Rfam, help='CM file with the Rfam library. Default %s''' % str(self.Rfam))
        ncRNA_group.add_argument('--cmsearch-CPUs', dest="cmsearch_CPUs", type=int, default = self.cmsearch_CPUs, help='''Number of CPUs to run cmsearch Default %s''' % str(self.cmsearch_CPUs))
        
####

    def check_parameters(self,args):
        """Check parameters consistency
            
        args -- set of parsed arguments"""

        if args.jsonFile==None:
            parser.print_help()
            sys.exit(-1)

        ##Getting outputs:
        working_dir = os.getcwd()

        if args.output_dir:
            self.output_dir = os.path.abspath(args.output_dir)
        else:
            self.output_dir = working_dir + "/step0" + str(args.annotation_step) + "_annotation_pipeline.V" + str(args.annotation_version)  + "/"

        if args.EVM_dir:
            self.EVM_dir = os.path.abspath(args.EVM_dir)
        else:
            self.EVM_dir = working_dir + "/step0" + str(int(args.annotation_step) + 1) + "_EVM.V" + str(args.annotation_version)  + "/"

        if args.dir_masked_chunks:
            self.dir_masked_chunks =os.path.abspath(args.dir_masked_chunks)
        else:
            self.dir_masked_chunks = self.output_dir + "/chunks_masked_reference"

        if args.dir_genome_chunks:
            self.dir_genome_chunks = os.path.abspath(args.dir_genome_chunks)
        else:
            self.dir_genome_chunks = self.output_dir + "/chunks_genome_reference"

        if args.dir_process_junctions:
            self.dir_process_junctions = os.path.abspath(args.dir_process_junctions)
        else:
            self.dir_process_junctions = self.output_dir + "/junctions"

        if args.augustus_prediction:
            self.augustus_prediction = os.path.abspath(args.augustus_prediction)
        else:
            self.augustus_prediction = self.output_dir + "/gene_predictions/augustus/augustus_gene_prediction.gff3"

        if args.augustus_preEVM:
            self.augustus_preEVM = os.path.abspath(args.augustus_preEVM)
        else:
            self.augustus_preEVM = self.output_dir + "/gene_predictions/augustus/augustus_preEVM.gff3" 

        if args.geneid_prediction:
            self.geneid_prediction = os.path.abspath(args.geneid_prediction)
        else: 
            self.geneid_prediction = self.output_dir + "/gene_predictions/geneid/geneid.gff3"

        if args.geneid_preEVM:
            self.geneid_preEVM = os.path.abspath(args.geneid_preEVM)
        else: 
            self.geneid_preEVM = self.output_dir + "/gene_predictions/geneid/geneid_preEVM.gff3"

        if args.genemark_prediction:
            self.genemark_prediction = os.path.abspath(args.genemark_prediction)
        else: 
            self.genemark_prediction = self.output_dir + "/gene_predictions/genemark.gtf"

        if args.genemark_preEVM:
            self.genemark_preEVM = os.path.abspath(args.genemark_preEVM)
        else: 
            self.genemark_preEVM = self.output_dir + "/gene_predictions/genemark_preEVM.gff3"

        if args.glimmer_prediction:
            self.glimmer_prediction = os.path.abspath(args.glimmer_prediction)
        else: 
            self.glimmer_prediction = self.output_dir + "/gene_predictions/glimmer.gff3"

        if args.glimmer_preEVM:
            self.glimmer_preEVM = os.path.abspath(args.glimmer_preEVM)
        else: 
            self.glimmer_preEVM = self.output_dir + "/gene_predictions/glimmer_preEVM.gff3"

        if args.spaln_cds:
            self.spaln_cds = os.path.abspath(args.spaln_cds)
        else: 
            self.spaln_cds =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_cds.gff3" 

        if args.spaln_gene:
            self.spaln_gene = os.path.abspath(args.spaln_gene)
        else: 
            self.spaln_gene =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_gene.gff3" 

        if args.geneid_introns_prediction:
            self.geneid_introns_prediction = os.path.abspath(args.geneid_introns_prediction)
        else: 
            self.geneid_introns_prediction = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns.gff3"

        if args.geneid_introns_preEVM:
            self.geneid_introns_preEVM = os.path.abspath(args.geneid_introns_preEVM)
        else: 
            self.geneid_introns_preEVM = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns_preEVM.gff3"

        if args.augustus_hints_prediction:
            self.augustus_hints_prediction = os.path.abspath(args.augustus_hints_prediction)
        else: 
            self.augustus_hints_prediction = self.output_dir + "/gene_predictions/augustus_with_hints/augustus_hints.gff3"

        if args.augustus_hints_preEVM:
            self.augustus_hints_preEVM = os.path.abspath(args.augustus_hints_preEVM)
        else: 
            self.augustus_hints_preEVM = self.output_dir + "/gene_predictions/augustus_with_hints/augustus_hints_preEVM.gff3"

        
        if args.genemark_ET_prediction:
            self.genemark_ET_prediction = os.path.abspath(args.genemark_ET_prediction)
        else: 
            self.genemark_ET_prediction = self.output_dir + "/gene_predictions/genemark-ET.gtf"

        if args.genemark_ET_preEVM:
            self.genemark_ET_preEVM = os.path.abspath(args.genemark_ET_preEVM)
        else: 
            self.genemark_ET_preEVM = self.output_dir + "/gene_predictions/genemark-ET_preEVM.gff3"

        if args.pasa_dir:
            self.pasa_dir = os.path.abspath(args.pasa_dir)
        else: 
            self.pasa_dir = self.output_dir + "/protein_and_transcript_mappings/pasa/" 

        if args.update_dir:
            self.update_dir = os.path.abspath(args.update_dir)
        else:
            self.update_dir = "step0" + str(int(args.annotation_step) + 2) + "_annotation_update.V" +str(args.annotation_version)  + "/" 

        if args.ncRNA_annotation_dir:
            self.ncRNA_annotation_dir = os.path.abspath(args.ncRNA_annotation_dir)
        else:
            self.ncRNA_annotation_dir = "step0" + str(int(args.annotation_step) + 3) + "_ncRNA_annotation.V" + str(args.annotation_version)  + "/"

        if args.out_cmsearch:
            self.out_cmsearch = os.path.abspath(args.out_cmsearch)
        else:
            self.out_cmsearch = self.ncRNA_annotation_dir + "/cmsearch.tbl"

        if args.out_tRNAscan:
            self.out_tRNAscan = os.path.abspath(args.out_tRNAscan)
        else:
            self.out_tRNAscan = self.ncRNA_annotation_dir + "/tRNAscan-SE/tRNAscan.out"

        ##Checking inputs
        if args.run_geneid or args.run_genemark or args.run_glimmer or args.run_augustus or args.run_geneid_introns or args.run_augustus_hints or args.run_genemark_ET:
            if args.genome_masked == None:
                print "Sorry! No masked genome fasta file defined"
            else:
                if not os.path.exists(args.genome_masked):
                    print args.genome_masked + " not found" 
                else:
                    args.genome_masked = os.path.abspath(args.genome_masked)

        if args.run_pasa or args.run_transdecoder or args.run_evm or args.run_update or args.run_ncRNA_annotation or args.run_spaln:
            if args.genome == None:
                print "Sorry! No genome fasta file defined"
            else:
                if not os.path.exists(args.genome):
                    print args.genome + " not found" 
                else:
                    args.genome = os.path.abspath(args.genome)

        if args.run_geneid_introns or args.run_augustus_hints or args.run_genemark_ET:
            if args.junctions == None and args.incoding_junctions == None:
                print "Sorry! No junctions gff file given."
            else:
                if args.junctions:
                    if not os.path.exists(args.junctions):
                        print args.junctions + " not found" 
                    else:
                        args.junctions = os.path.abspath(args.junctions)
                if args.incoding_junctions:
                    if os.path.exists(args.incoding_junctions):
                        args.incoding_junctions = os.path.abspath(args.incoding_junctions)      

        if args.run_augustus or args.run_augustus_hints:
            if args.species == None:
                print "Sorry! No species for augustus defined"        

        if args.run_geneid or args.run_update or args.run_geneid_introns:
            if args.geneid_parameters == None:
                print "Sorry! No geneid parameters file defined" 
            else:
                if not os.path.exists(args.geneid_parameters):
                    print args.geneid_parameters + " not found" 
                else:
                    args.geneid_parameters = os.path.abspath(args.geneid_parameters) 

        if args.run_glimmer:
            if args.glimmer_directory == None:
                print "Sorry! No glimmer trained directory given"
            else:
                if not os.path.exists(args.glimmer_directory):
                    print args.glimmer_directory + " not found" 
                else:
                    args.glimmer_directory = os.path.abspath(args.glimmer_directory)

        if args.run_spaln or args.run_ncRNA_annotation:
            if args.proteins == None:
                print "Sorry! No protein evidence file given."
            else:
                if not os.path.exists(args.proteins):
                    print args.proteins + " not found" 
                else:
                    args.proteins = os.path.abspath(args.proteins)

        if args.run_augustus_hints:
            if args.extrinsic_file_augustus_hints == None:
                print "Sorry! No extrinsic file for augustus with hints given."
            else:
                if not os.path.exists(args.extrinsic_file_augustus_hints):
                    print args.extrinsic_file_augustus_hints + " not found" 
                else:
                    args.extrinsic_file_augustus_hints = os.path.abspath(args.extrinsic_file_augustus_hints)

            if args.ep_hints:
                if not os.path.exists(args.ep_hints):
                    print args.ep_hints + "not found"
                else:
                    args.ep_hints = os.path.abspath(args.ep_hints) 

        if args.run_pasa or args.run_transdecoder or args.run_update:
            if args.pasadb == None:
                print "Sorry! No pasadb name given"
            
            if not os.path.exists(args.pasa_home):
                print args.pasa_home + " not found" 
            else:
                self.pasa_home = os.path.abspath(args.pasa_home)

        if args.run_pasa:
            if args.transcripts == None:
                print "Sorry! No transcript evidence file found."
            else:
                if not os.path.exists(args.transcripts):
                    print args.transcripts + " not found" 
                else:
                    args.transcripts = os.path.abspath(args.transcripts) 
            
            if args.pasa_config == None:
                print "Sorry! No pasa configuration file found."
            else:
                if not os.path.exists(args.pasa_config):
                    print args.pasa_config + " not found" 
                else:
                    args.pasa_config = os.path.abspath(args.pasa_config)

            if args.create_database:
                print "WARNING: Remember to turn the create-database parameter to false if it's not the first time that you run PASA!"
            else:
                print "WARNING: Remember to set the create-database parameter to true the first time you run PASA!"

        if args.cufflinks:
            if not os.path.exists(args.cufflinks):
                print args.cufflinks + " not found" 
            else:
                args.cufflinks = os.path.abspath(args.cufflinks)

        if not os.path.exists(args.evm_script):
            print args.evm_script + " not found" 
        else:
            args.evm_script = os.path.abspath(args.evm_script)

        if args.run_update:
            if args.update_config == None:
                print "Sorry! No pasa update configuration file found."
            else:
                if not os.path.exists(args.update_config):
                    print args.update_config + " not found" 
                else:
                    args.update_config = os.path.abspath(args.update_config)

        if not os.path.exists(args.Rfam):
            print args.Rfam + " not found " 

        if args.run_ncRNA_annotation:
            if not args.ncRNA_version:
                print "Sorry! ncRNA annotation version given!"

            if args.RM_gff == None:
                print "Sorry! No repeat masker gff output file found."
            else:
                if not os.path.exists(args.RM_gff):
                    print args.RM_gff + " not found" 
                else:
                    args.RM_gff = os.path.abspath(args.RM_gff)
        

        if args.run_update or args.run_ncRNA_annotation:
            if args.project_name == None:
                print "Sorry! No project name and annotation version found."
   
###

    def storeGeneralParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.generalParameters["jsonFile"] = args.jsonFile
        self.generalParameters["run_geneid"] = args.run_geneid
        self.generalParameters["run_augustus"] = args.run_augustus
        self.generalParameters["run_genemark"] = args.run_genemark
        self.generalParameters["run_glimmer"] = args.run_glimmer
        self.generalParameters["run_geneid_introns"] = args.run_geneid_introns
        self.generalParameters["run_augustus_hints"] = args.run_augustus_hints
        self.generalParameters["run_genemark_ET"] = args.run_genemark_ET
        self.generalParameters["run_spaln"] = args.run_spaln
        self.generalParameters["run_pasa"] = args.run_pasa
        self.generalParameters["run_transdecoder"] = args.run_transdecoder
        self.generalParameters["run_evm"] = args.run_evm
        self.generalParameters["run_update"] = args.run_update
        self.generalParameters["run_ncRNA_annotation"] = args.run_ncRNA_annotation
        self.generalParameters["pipeline_HOME"] = args.pipeline_HOME
        self.allParameters  ["Parameters"] = self.generalParameters

    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.inputParameters["genome"] = args.genome
        self.inputParameters["genome_masked"] = args.genome_masked
        self.inputParameters["junctions"] = args.junctions
        self.inputParameters["incoding_junctions"] = args.incoding_junctions
        self.inputParameters["species"] = args.species
        self.inputParameters["geneid_parameters"] = args.geneid_parameters
        self.inputParameters["glimmer_directory"] = args.glimmer_directory
        self.inputParameters["proteins"] = args.proteins
        self.inputParameters["extrinsic_file_augustus_hints"] = args.extrinsic_file_augustus_hints
        self.inputParameters["ep_hints"] = args.ep_hints
        self.inputParameters["pasadb"] = args.pasadb
        self.inputParameters["transcripts"] = args.transcripts
        self.inputParameters["pasa_config"] = args.pasa_config
        self.inputParameters["RM_gff"] = args.RM_gff
        self.inputParameters["cufflinks"] = args.cufflinks
        self.inputParameters["update_config"] = args.update_config
        self.inputParameters["project_name"] = args.project_name
        self.inputParameters["ncRNA_version"] = args.ncRNA_version
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates output parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """

        self.outputParameters["annotation_step"] = args.annotation_step
        self.outputParameters["annotation_version"] = args.annotation_version
        self.outputParameters["output_dir"] = self.output_dir
        self.outputParameters ["dir_masked_chunks"] = self.dir_masked_chunks
        self.outputParameters ["dir_genome_chunks"] = self.dir_genome_chunks
        self.outputParameters ["dir_process_junctions"] = self.dir_process_junctions
        self.outputParameters["EVM_dir"] = self.EVM_dir
        self.outputParameters["augustus_prediction"] = self.augustus_prediction
        self.outputParameters["augustus_preEVM"] = self.augustus_preEVM
        self.outputParameters["geneid_prediction"] = self.geneid_prediction
        self.outputParameters["geneid_preEVM"] = self.geneid_preEVM
        self.outputParameters["genemark_prediction"] = self.genemark_prediction
        self.outputParameters["genemark_preEVM"] = self.genemark_preEVM
        self.outputParameters["glimmer_prediction"] = self.glimmer_prediction
        self.outputParameters["glimmer_preEVM"] = self.glimmer_preEVM
        self.outputParameters["spaln_cds"] = self.spaln_cds
        self.outputParameters["spaln_gene"] = self.spaln_gene
        self.outputParameters["geneid_introns_prediction"] = self.geneid_introns_prediction
        self.outputParameters["geneid_introns_preEVM"] = self.geneid_introns_preEVM
        self.outputParameters["augustus_hints_prediction"] = self.augustus_hints_prediction
        self.outputParameters["augustus_hints_preEVM"] = self.augustus_hints_preEVM
        self.outputParameters["genemark_ET_prediction"] = self.genemark_ET_prediction
        self.outputParameters["genemark_ET_preEVM"] = self.genemark_ET_preEVM
        self.outputParameters["pasa_dir"] = self.pasa_dir
        self.outputParameters["update_dir"] = self.update_dir
        self.outputParameters["ncRNA_annotation_dir"] = self.ncRNA_annotation_dir
        self.outputParameters["out_cmsearch"] = self.out_cmsearch
        self.outputParameters["out_tRNAscan"] = self.out_tRNAscan
        self.allParameters["Outputs"] = self.outputParameters

    def storeChunksParameters(self,args):
        """Updates chunks parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """

        self.chunksParameters["masked_chunks"] = args.masked_chunks
        self.chunksParameters["genome_chunks"] = args.genome_chunks
        self.chunksParameters["protein_chunks"] = args.protein_chunks
        self.allParameters["Chunks"] = self.chunksParameters
     
    def storeAugustusParameters(self,args):
        """Updates augustus parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusParameters["aug_alternatives_from_sampling"] = args.aug_alternatives_from_sampling
        self.augustusParameters["aug_uniqueGeneId"] = args.aug_uniqueGeneId
        self.augustusParameters["aug_gff3"] = args.aug_gff3
        self.augustusParameters["aug_sample"] = args.aug_sample
        self.augustusParameters["aug_noInFrameStop"] = args.aug_noInFrameStop
        self.augustusParameters["aug_maxtracks"] = args.aug_maxtracks
        self.augustusParameters["aug_singlestrand"] = args.aug_singlestrand
        self.augustusParameters["aug_strand"] = args.aug_strand
        self.augustusParameters["aug_min_intron_len"] = args.aug_min_intron_len
        self.augustusParameters["augustus_weights"] = args.augustus_weights
        self.augustusParameters["additional_augustus_options"] = args.additional_augustus_options    
        self.allParameters ["augustus"] = self.augustusParameters

    def storeGeneidParameters(self,args):
        """Updates geneid parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.geneidParameters["geneid_weights"] = args.geneid_weights
        self.geneidParameters["geneid_options"] = args.geneid_options       
        self.allParameters ["geneid"] = self.geneidParameters
        
    def storeGenemarkParameters(self,args):
        """Updates genemark parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genemarkParameters["gmk_max_contig"] = args.gmk_max_contig
        self.genemarkParameters["gmk_min_contig"] = args.gmk_min_contig
        self.genemarkParameters["gmk_max_gap"] = args.gmk_max_gap
        self.genemarkParameters["gmk_cores"] = args.gmk_cores
        self.genemarkParameters["genemark_weights"] = args.genemark_weights
        self.genemarkParameters["additional_genemark_options"] = args.additional_genemark_options   
        self.allParameters ["genemark"] = self.genemarkParameters

    def storeGlimmerParameters(self,args):
        """Updates glimmer parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.glimmerParameters["glimmer_weights"] = args.glimmer_weights
        self.allParameters ["glimmer"] = self.glimmerParameters

    def storeSpalnParameters(self,args):
        """Updates spaln parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.spalnParameters["spaln_ya"] = args.spaln_ya
        self.spalnParameters["spaln_M"] = args.spaln_M
        self.spalnParameters["spaln_O"] = args.spaln_O
        self.spalnParameters["spaln_Q"] = args.spaln_Q
        self.spalnParameters["spaln_t"] = args.spaln_t
        self.spalnParameters["spaln_weights"] = args.spaln_weights
        self.spalnParameters["additional_spaln_options"] = args.additional_spaln_options      
        self.allParameters ["spaln"] = self.spalnParameters

    def storeGeneidIntronsParameters(self,args):
        """Updates geneid with introns parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.geneidIntronsParameters["geneid_introns_weights"] = args.geneid_introns_weights
        self.geneidIntronsParameters["geneid_introns_options"] = args.geneid_introns_options       
        self.allParameters ["geneid_introns"] = self.geneidIntronsParameters

    def storeAugustusHintsParameters(self,args):
        """Updates augustus with introns parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusHintsParameters["augustus_hints_weights"] = args.augustus_hints_weights
        self.augustusHintsParameters["additional_augustus_hints_options"] = args.additional_augustus_hints_options       
        self.allParameters ["augustus_hints"] = self.augustusHintsParameters

    def storeGenemarkETParameters(self,args):
        """Updates genemark parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genemarkETParameters["genemark_ET_weights"] = args.genemark_ET_weights
        self.genemarkETParameters["gmk_et_score"] = args.gmk_et_score
        self.genemarkETParameters["additional_genemark_ET_options"] = args.additional_genemark_ET_options 
        self.allParameters ["genemark-ET"] = self.genemarkETParameters

    def storePasaParameters(self,args):
        """Updates pasa parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pasaParameters["pasa_CPU"] = args.pasa_CPU
        self.pasaParameters["create_database"] = args.create_database
        self.pasaParameters["pasa_step"] = args.pasa_step
        self.pasaParameters["pasa_weights"] = args.pasa_weights
        self.pasaParameters["pasa_home"] = args.pasa_home
        self.allParameters["pasa"] = self.pasaParameters

    def storeTransdecoderParameters(self,args):
        """Updates transdecoder parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.transdecoderParameters["transdecoder_weights"] = args.transdecoder_weights
        self.allParameters["transdecoder"] = self.transdecoderParameters

    def storeEvmParameters(self,args):
        """Updates evm parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.evmParameters["evm_script"] = args.evm_script
        self.allParameters["evm"] = self.evmParameters

    def storencRNAannotationParameters(self,args):
        """Updates ncRNA Annotation parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.ncRNAannotationParameters["cmsearch_CPUs"] = args.cmsearch_CPUs
        self.ncRNAannotationParameters["Rfam"] = args.Rfam
        self.allParameters["ncRNA_annotation"] = self.ncRNAannotationParameters

#####

#1.Create object class Configuration File
configManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="create_configuration_file",
                description="Create a configuration json file for the annotation pipeline."
                )     

#2.1 Updates arguments and parsing
configManager.register_parameter(parser)

args = parser.parse_args()

#2.2 Check Parameters
configManager.check_parameters(args)

#3. store arguments to super map structure
configManager.storeGeneralParameters(args)
configManager.storeInputParameters(args)
configManager.storeOutputParameters(args)
configManager.storeChunksParameters(args)
configManager.storeAugustusParameters(args)
configManager.storeGeneidParameters(args)
configManager.storeGenemarkParameters(args)
configManager.storeGlimmerParameters(args)
configManager.storeSpalnParameters(args)
configManager.storeGeneidIntronsParameters(args)
configManager.storeAugustusHintsParameters(args)
configManager.storeGenemarkETParameters(args)
configManager.storePasaParameters(args)
configManager.storeTransdecoderParameters(args)
configManager.storeEvmParameters(args)
configManager.storencRNAannotationParameters(args)

###

#4. Store JSON file
with open(args.jsonFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)

