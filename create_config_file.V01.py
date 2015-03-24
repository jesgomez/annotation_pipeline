#!/usr/bin/env python
import os
import json
import argparse
import sys

class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor"""
        #GENERAL PARAMETERS
        self.run_geneid = True                             #By default run geneid step
        self.run_augustus = True                           #By default run augustus step
        self.run_genemark = True                           #By default run genemark step
        self.run_glimmer = True                            #By default run glimmer step
        self.run_geneid_introns = True                     #By default run geneid with introns step
        self.run_augustus_introns = True                   #By default run augustus with introns step
        self.run_genemark_ET = True                        #By default run genemark-ET step
        self.run_spaln = True                              #By default run spaln step
        self.run_pasa = True                               #By default run pasa step
        self.run_transdecoder = True                       #By default run transdecoder step
        self.run_evm = True                                #By default run EvidenceModeler step
        self.run_update = True                             #By default run annotation update step
        self.pipeline_HOME = "/project/devel/aateam/src/Annotation_pipeline.V01/"   # Path to the pipeline home directory
 
        #INPUT PARAMETERS
        self.jsonFile = None                               #Name of the json configuration file to be created.
        self.genome = None                                 #Path to the fasta genome.
        self.genome_masked = None                          #Path to the fasta genome masked.
        self.species = None                                #Species name to run augustus with its trained parameters.
        self.geneid_parameters = None                      #Path to the geneid parameters file.
        self.glimmer_directory = None                      #Path to the directory containing the trained parameters for running glimer.
        self.junctions = None                              #Path to the junctions gff file to run gene predictors with introns.
        self.incoding_junctions = None                     #Path to the junctions in coding regions gff file to run gene predictors with introns.
        self.extrinsic_file_augustus_introns = None        #Extrinsic file to use when running augustus with introns. For more information read augustus documentation. 
        self.proteins = None                               #Path to the fasta with protein evidence.
        self.pasadb = None                                 #Name of the pasa database, it must coincide with the name given in pasa-config.
        self.transcripts = None                            #Path to the fasta with transcript evidence.
        self.pasa_config = None                            #Path to the Pasa configuration file.
        self.cufflinks = None                              #Path to the cufflinks gtf file.
        self.update_config = None                          #Path to the Pasa configuration file for the annotation update step.
        self.project_name = None                           #Name of the project and version, to give the names to the final annotation output. 

        #OUTPUT PARAMETERS
        self.annotation_step = "3"                        #Step of the annotation pipeline in the annotation process.    
        self.annotation_version = "01"                     #Version of the annotation process.      
        self.output_dir = "step0" + self.annotation_step + "_annotation_pipeline.V" + self.annotation_version  + "/"   #Directory to keep the outputs of the first annotation steps.   
        self.EVM_inputs = "step0" + str(int(self.annotation_step) + 1) + "_EVM.V" + self.annotation_version  + "/"   #Directory to keep the files for the EVM step.    
        self.update_dir = "step0" + str(int(self.annotation_step) + 2) + "_annotation_update.V" + self.annotation_version  + "/"   #Directory to keep the files for annotation update step.     
        self.augustus_prediction = self.output_dir + "/gene_predictions/augustus/augustus_gene_prediction.gff3"  #Output file for the augustus predictions.
        self.augustus_preEVM = self.output_dir + "/gene_predictions/augustus/augustus_preEVM.gff3"  #Output file for the augustus predictions converted for EVM.
        self.geneid_prediction = self.output_dir + "/gene_predictions/geneid/geneid.gff3"  #Output file for the geneid predictions.
        self.geneid_preEVM = self.output_dir + "/gene_predictions/geneid/geneid_preEVM.gff3"  #Output file for the geneid predictions converted for EVM.
        self.genemark_prediction = self.output_dir + "/gene_predictions/genemark.gtf"  #Output file for the genemark predictions.
        self.genemark_preEVM = self.output_dir + "/gene_predictions/genemark_preEVM.gff3"  #Output file for the genemark predictions converted for EVM.
        self.glimmer_prediction = self.output_dir + "/gene_predictions/glimmer.gff3"  #Output file for the glimmer predictions.
        self.glimmer_preEVM = self.output_dir + "/gene_predictions/glimmer_preEVM.gff3"  #Output file for the glimmer predictions converted for EVM.
        self.geneid_introns_prediction = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns.gff3"  #Output file for the geneid with introns predictions.
        self.geneid_introns_preEVM = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns_preEVM.gff3"  #Output file for the geneid with introns predictions converted for EVM.
        self.augustus_introns_prediction = self.output_dir + "/gene_predictions/augustus_with_introns/augustus_introns.gff3"  #Output file for the augustus with introns predictions.
        self.augustus_introns_preEVM = self.output_dir + "/gene_predictions/augustus_with_introns/augustus_introns_preEVM.gff3"  #Output file for the augustus with introns predictions converted for EVM.
        self.genemark_ET_prediction = self.output_dir + "/gene_predictions/genemark-ET.gtf"  #Output file for the genemark-ET predictions.
        self.genemark_ET_preEVM = self.output_dir + "/gene_predictions/genemark-ET_preEVM.gff3"  #Output file for the genemark-ET predictions converted for EVM.
        self.spaln_gene =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_gene.gff3"    #Output file for the spaln output in a gene gff3 format.
        self.spaln_cds =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_cds.gff3"    #Output file for the spaln output in a cds gff3 format.
        self.pasa_dir =  self.output_dir + "/protein_and_transcript_mappings/pasa/"     #Directory to keep all the pasa outputs.         

        #AUGUSTUS PARAMETERS
        self.alternatives_from_sampling = "true"           #Report alternative transcripts generated through probabilistic sampling.                     
        self.uniqueGeneId = "true"                         #If true, output gene identifyers like this: seqname.gN.         
        self.gff3 = "ON"                                   #Output in gff3 format.
        self.sample = 60 
        self.noInFrameStop = "true"                        #Do not report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur.
        self.maxtracks = 2                                 #Maximum number of tracks allowed. 
        self.singlestrand = "false"                        #Predict genes independently on each strand, allow overlapping genes on opposite strands.
        self.strand= "both"    
        self.min_intron_len= 30                            #Minimum predicted intron length.                            
        self.augustus_weights = [2, 2, 1]                  #Weights given to augustus predictions when running EVM.
        self.additional_augustus_options = None            #Additional augustus options to run it, see augustus help for more information.

        #GENEID PARAMETERS
        self.geneid_weights = [2,1,2]                      #Weights given to geneid predictions when running EVM.
        self.geneid_options = "3U "                      #Desired geneid options to run it, see geneid documentation for more information.
  
        #GENEMARK PARAMETERS
        self.min_contig = 50000                            #Will ignore contigs shorter then min_contig in training
        self.max_contig = 5000000                          #will split input genomic sequence into contigs shorter than max_contig.
        self.max_gap = 5000                                #Will split sequence at gaps longer than max_gap. Letters 'n' and 'N' are       interpreted as standing within gaps 
        self.genemark_cores = 8                                     #Number of threads for running genemark. 
        self.additional_genemark_options = None            #Additional genemark options to run it, see genemark documentation for more information.
        self.genemark_weights = [1, 1, 1]                  #Weights given to genemark predictions when running EVM.

        #GLIMMER PARAMETERS
        self.glimmer_weights = [1, 1, 1]                   #Weights given to glimmer predictions when running EVM.

        #GENEID INTRONS PARAMETERS
        self.geneid_introns_weights = [3,3,3]              #Weights given to geneid with intron predictions when running EVM.
        self.geneid_introns_options = "3nU"                #Desired geneid options to run geneid with introns, see geneid documentation for more information.

        #AUGUSTUS INTRONS PARAMETERS
        self.augustus_introns_weights = [3,3,3]              #Weights given to augustus with intron predictions when running EVM.
        self.additional_augustus_introns_options = None                 #Desired augustus options to run geneid with introns, see geneid documentation for more information.

        #GENEMARK-ET PARAMETERS
        self.et_score = 4                                 #Minimum score of intron in initiation of the ET algorithm
        self.genemark_ET_weights = [3,3,3]                #Weights given to augustus with intron predictions when running EVM.
        self.additional_genemark_ET_options = None        #Additional genemark-ET options to run it, see genemark documentation for more information.

        #SPALN PARAMETERS
        self.ya = 1                                        #Stringency of splice site. 0->3: strong->weak.
        self.M = 4                                         #Number of outputs per query (1) (max=4 for genome vs cDNA|protein).
        self.O = 0                                         #0:Gff3_gene; 1:alignment; 2:Gff3_match; 3:Bed; 4:exon-inf; 5:intron-inf; 6:cDNA; 7:translated; 8: block-only; 12: binary.
        self.Q = 7                                         # 0:DP; 1-3: HSP-Search; 4-7; Block-Search.
        self.t = 8                                         #Number of threads.
        self.spaln_weights = [10, 8, 10]                   #Weights given to spaln mappings when running EVM.
        self.additional_spaln_options = None               #Additional spaln options to run it, see spaln help for more information.

        #PASA PARAMETERS
        self.pasa_CPU = 8                                       #Number of pasa_CPUs to run Pasa.
        self.pasa_step = 1                                 #Step from where to start running Pasa.
        self.pasa_weights = [8, 10, 8]                     #Weights given to pasa mappings when running EVM. 
        self.create_database = False                       #By default do not create pasa database.
  
        #TRANSDECODER PARAMETERS
        self.transdecoder_weights = [3, 2, 3]              #Weights given to the pasa transdecoder output gff3 file.

        #EVM PARAMETERS
        self.evm_script = "/project/devel/aateam/bin/annotation_scripts/evm.sh"

        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.augustusParameters = {}
        self.geneidParameters = {}
        self.genemarkParameters = {}
        self.glimmerParameters = {}
        self.geneidIntronsParameters = {}
        self.augustusIntronsParameters = {}
        self.genemarkETParameters = {}
        self.spalnParameters = {}
        self.pasaParameters = {}
        self.transdecoderParameters = {}
        self.evmParameters = {}

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_augustus(parser)
        self.register_geneid(parser)
        self.register_genemark(parser)
        self.register_glimmer(parser)
        self.register_geneid_introns(parser)
        self.register_augustus_introns(parser)
        self.register_genemark_ET(parser)
        self.register_spaln(parser)
        self.register_pasa(parser)
        self.register_transdecoder(parser)
        self.register_evm(parser)

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('--no-geneid', dest="run_geneid", action="store_false", help='If specified, do not run geneid step.')
        general_group.add_argument('--no-augustus', dest="run_augustus", action="store_false", help='If specified, do not run augustus step.')
        general_group.add_argument('--no-genemark', dest="run_genemark", action="store_false", help='If specified, do not run genemark step.')
        general_group.add_argument('--no-glimmer', dest="run_glimmer", action="store_false", help='If specified, do not run glimmer step.')
        general_group.add_argument('--no-geneid-introns', dest="run_geneid_introns", action="store_false", help='If specified, do not run geneid with introns step.')
        general_group.add_argument('--no-augustus-introns', dest="run_augustus_introns", action="store_false", help='If specified, do not run augustus with introns step.')
        general_group.add_argument('--no-genemark-ET', dest="run_genemark_ET", action="store_false", help='If specified, do not run genemark-ET step.')
        general_group.add_argument('--no-spaln', dest="run_spaln", action="store_false", help='If specified, do not run spaln step.')
        general_group.add_argument('--no-pasa', dest="run_pasa", action="store_false", help='If specified, do not run pasa step.')
        general_group.add_argument('--no-transdecoder', dest="run_transdecoder", action="store_false", help='If specified, do not run transdecoder step.')
        general_group.add_argument('--no-evm', dest="run_evm", action="store_false", help='If specified, do not run EVM step.')
        general_group.add_argument('--no-update', dest="run_update", action="store_false", help='If specified, do not run the annotation update step.')
        general_group.add_argument('--pipeline-HOME', dest="pipeline_HOME", metavar="pipeline_HOME", help='Path to the pipeline home directory. Default %s' %self.pipeline_HOME)

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--genome', dest="genome", metavar="genome", help='Path to the fasta genome.')
        input_group.add_argument('--genome-masked', dest="genome_masked", metavar="genome_masked", help='Path to the fasta genome masked.')
        input_group.add_argument('-jsonFile', dest="jsonFile", metavar="jsonFile", help='Configuration JSON to be generated. Default %s' % self.jsonFile)
        input_group.add_argument('--species', dest="species", metavar="species", help='Species name to run augustus with its trained parameters.')
        input_group.add_argument('--geneid-parameters', dest="geneid_parameters", metavar="geneid_parameters", help='Path to the geneid parameters file.')
        input_group.add_argument('--glimmer-directory', dest="glimmer_directory", metavar="glimmer_directory", help='Path to the directory containing the trained parameters for running glimer.')
        input_group.add_argument('--junctions', dest="junctions", metavar="junctions", help='Path to the junctions gff file to run gene predictors with introns. Do not needed if incoding junctions existent file is given')
        input_group.add_argument('--incoding-junctions', dest="incoding_junctions", metavar="incoding_junctions", help='Path to the junctions in coding regions gff file to run gene predictors with introns. (Optional, it can be given if this file already exists or if you want to keep the processed incoding junctions in that concrete path.)')
        input_group.add_argument('--extrinsic-file-augustus-introns', dest="extrinsic_file_augustus_introns", metavar="extrinsic_file_augustus_introns", help='Path to the Extrinsic file to use when running augustus with introns. For more information read augustus documentation.')
        input_group.add_argument('--proteins', dest="proteins", metavar="proteins", help='Path to the fasta with protein evidence.')
        input_group.add_argument('--pasadb', dest="pasadb", metavar="pasadb", help='Name of the pasa database, it must coincide with the name in pasa_config.')
        input_group.add_argument('--transcripts', dest="transcripts", metavar="transcripts", help='Path to the fasta with transcript evidence.')
        input_group.add_argument('--pasa-config', dest="pasa_config", metavar="pasa_config", help='Path to the pasa configuration file.')
        input_group.add_argument('--cufflinks', dest="cufflinks", metavar="cufflinks", help='Path to the cufflinks gtf file.')
        input_group.add_argument('--update-config', dest="update_config", metavar="update_config", help='Path to the Pasa configuration file.')
        input_group.add_argument('--project-name', dest="project_name", metavar="project_name", nargs="+", help='Name of the project and version of the annotation space separated, to give the names to the final annotation output.')

    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--annotation-step', dest="annotation_step", type=int, help='Step of the annotation pipeline in the annotation process. Default %s' % self.annotation_step)
        output_group.add_argument('--annotation-version', dest="annotation_version", help='Version of the annotation process. Default %s' % self.annotation_version)
        output_group.add_argument('--output-dir', dest="output_dir", help='Directory to keep the outputs of the first annotation steps.  Default %s' % self.output_dir)
        output_group.add_argument('--EVM-inputs', dest="EVM_inputs", help='Directory to keep the files for the EVM step. Default %s' % self.EVM_inputs)
        output_group.add_argument('--augustus-prediction', dest="augustus_prediction", help='Output file for the augustus predictions.  Default %s' % self.augustus_prediction)
        output_group.add_argument('--augustus-preEVM', dest="augustus_preEVM", help='Output file for the augustus predictions converted for EVM.  Default %s' % self.augustus_preEVM)
        output_group.add_argument('--geneid-prediction', dest="geneid_prediction", help='Output file for the geneid predictions.  Default %s' % self.geneid_prediction)
        output_group.add_argument('--geneid-preEVM', dest="geneid_preEVM", help='Output file for the geneid predictions converted for EVM.  Default %s' % self.geneid_preEVM)
        output_group.add_argument('--genemark-prediction', dest="genemark_prediction", help='Output file for the genemark predictions.  Default %s' % self.genemark_prediction)
        output_group.add_argument('--genemark-preEVM', dest="genemark_preEVM", help='Output file for the genemark predictions converted for EVM.  Default %s' % self.genemark_preEVM)
        output_group.add_argument('--glimmer-prediction', dest="glimmer_prediction", help='Output file for the glimmer predictions.  Default %s' % self.glimmer_prediction)
        output_group.add_argument('--glimmer-preEVM', dest="glimmer_preEVM", help='Output file for the glimmer predictions converted for EVM.  Default %s' % self.glimmer_preEVM)
        output_group.add_argument('--geneid-introns-prediction', dest="geneid_introns_prediction", help='Output file for the geneid with introns predictions.  Default %s' % self.geneid_introns_prediction)
        output_group.add_argument('--geneid-introns-preEVM', dest="geneid_introns_preEVM", help='Output file for the geneid with introns predictions converted for EVM.  Default %s' % self.geneid_preEVM)
        output_group.add_argument('--augustus-introns-prediction', dest="augustus_introns_prediction", help='Output file for the augustus with introns predictions.  Default %s' % self.augustus_introns_prediction)
        output_group.add_argument('--augustus-introns-preEVM', dest="augustus_introns_preEVM", help='Output file for the augustus with introns predictions converted for EVM.  Default %s' % self.augustus_preEVM)
        output_group.add_argument('--genemark-ET-prediction', dest="genemark_ET_prediction", help='Output file for the genemark-ET predictions.  Default %s' % self.genemark_ET_prediction)
        output_group.add_argument('--genemark-ET-preEVM', dest="genemark_ET_preEVM", help='Output file for the genemark-ET predictions converted for EVM.  Default %s' % self.genemark_ET_preEVM)
        output_group.add_argument('--spaln-gene', dest="spaln_gene", help='Output file for the spaln output in a gene gff3 format.  Default %s' % self.spaln_gene)
        output_group.add_argument('--spaln-cds', dest="spaln_cds", help='Output file for the spaln output in a cds gff3 format.  Default %s' % self.spaln_cds)
        output_group.add_argument('--pasa-dir', dest="pasa_dir", help='Directory to keep all the pasa outputs.  Default %s' % self.pasa_dir)
        output_group.add_argument('--update-dir', dest="update_dir", help='Directory to keep the files for annotation update step.   Default %s' % self.update_dir)

    def register_augustus(self, parser):
        """Register all augustus parameters with given
        argparse parser

        parser -- the argparse parser
        """
        augustus_group = parser.add_argument_group('Augustus parameters')
        augustus_group.add_argument('--alternatives-from-sampling', dest="alternatives_from_sampling", help='''Report alternative transcripts generated through probabilistic, possible values true/false. Default %s''' % str(self.alternatives_from_sampling))
        augustus_group.add_argument('--uniqueGeneId', dest="uniqueGeneId", help='''If true, output gene identifyers like this: seqname.gN, possible values true/false. Default %s''' % str(self.uniqueGeneId))
        augustus_group.add_argument('--gff3', dest="gff3", help='''Output in gff3 format, possible values ON/OFF. Default %s''' % str(self.gff3))
        augustus_group.add_argument('--sample', dest="sample", type=int, help='''Default %s''' % str(self.sample))
        augustus_group.add_argument('--noInFrameStop', dest="noInFrameStop", help='''Do not report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur, possible values true/false. Default %s''' % str(self.noInFrameStop))
        augustus_group.add_argument('--maxtracks', dest="maxtracks", type=int, help='''Maximum number of tracks allowed. Default %s''' % str(self.maxtracks))
        augustus_group.add_argument('--singlestrand', dest="singlestrand", help='''Predict genes independently on each strand, allow overlapping genes on opposite strands. possible values true/false. Default %s''' % str(self.singlestrand))
        augustus_group.add_argument('--strand', dest="strand", help='''Possible values both/forward/backward. Default %s''' % str(self.strand))
        augustus_group.add_argument('--min-intron-len', dest="min_intron_len", type=int, help='''Default %s''' % str(self.min_intron_len))
        augustus_group.add_argument('--augustus-weights', dest="augustus_weights", nargs="+", help='Weights given to augustus predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 2 2 1')
        augustus_group.add_argument('--additional-augustus-options', dest="additional_augustus_options", help='Additional augustus options to run it, see augustus help for more information about the possible options.')
       
    def register_geneid(self, parser):
        """Register all geneid parameters with given
        argparse parser

        parser -- the argparse parser
        """
        geneid_group = parser.add_argument_group('Geneid parameters')
        geneid_group.add_argument('--geneid-weights', dest="geneid_weights", nargs="+", help='Weights given to geneid predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 2 1 2 ')
        geneid_group.add_argument('--geneid-options', dest="geneid_options", help='Desired geneid options to run it, see geneid documentation for more information. Default %s''' % str(self.geneid_options))

    def register_genemark(self, parser):
        """Register all genemark parameters with given
        argparse parser

        parser -- the argparse parser
        """
        genemark_group = parser.add_argument_group('Genemark parameters')
        genemark_group.add_argument('--min-contig', dest="min_contig", type=int, help='''Will ignore contigs shorter then min_contig in training. Default %s''' % str(self.min_contig))
        genemark_group.add_argument('--max-contig', dest="max_contig", type=int, help='''Will split input genomic sequence into contigs shorter than max_contig. Default %s''' % str(self.max_contig))
        genemark_group.add_argument('--max-gap', dest="max_gap", type=int, help='''Will split sequence at gaps longer than max_gap. Letters 'n' and 'N' are interpreted as standing within gaps. Default %s''' % str(self.max_gap))
        genemark_group.add_argument('--genemark-cores', dest="genemark_cores", type=int, help='''Number of threads for running genemark. Default %s''' % str(self.genemark_cores))
        genemark_group.add_argument('--additional-genemark-options', dest="additional_genemark_options", help='Additional genemark options to run it, see genemark documentation for more information.')
        genemark_group.add_argument('--genemark-weights', dest="genemark_weights", nargs="+", help='Weights given to genemark predictions when running EVM. Specify the weight for each EVM run separated by a space. Example  1 1 1 ') 

    def register_glimmer(self, parser):
        """Register all glimmer parameters with given
        argparse parser

        parser -- the argparse parser
        """
        glimmer_group = parser.add_argument_group('Glimmer parameters')
        glimmer_group.add_argument('--glimmer-weights', dest="glimmer_weights", nargs="+", help='Weights given to glimmer predictions when running EVM. Specify the weight for each EVM run separated by a space. Example  1 1 1 ')

    def register_geneid_introns(self, parser):
        """Register all geneid with introns parameters with given
        argparse parser

        parser -- the argparse parser
        """
        geneid_introns_group = parser.add_argument_group('Geneid Introns parameters')
        geneid_introns_group.add_argument('--geneid-introns-weights', dest="geneid_introns_weights", nargs="+", help='Weights given to geneid with intron predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 3 3 3 ')
        geneid_introns_group.add_argument('--geneid-introns-options', dest="geneid_introns_options", help='Desired geneid with intron options to run it, see geneid documentation for more information. Default %s''' % str(self.geneid_options))

    def register_augustus_introns(self, parser):
        """Register all augustus with introns parameters with given
        argparse parser

        parser -- the argparse parser
        """
        augustus_introns_group = parser.add_argument_group('Augustus Introns parameters')
        augustus_introns_group.add_argument('--augustus-introns-weights', dest="augustus_introns_weights", nargs="+", help='Weights given to augustus with intron predictions when running EVM. Specify the weight for each EVM run separated by a space. Example 3 3 3 ')
        augustus_introns_group.add_argument('--additional-augustus-introns-options', dest="additional_augustus_introns_options", help='Desired augustus with intron options to run it, see augustus documentation for more information.''')

    def register_genemark_ET(self, parser):
        """Register all genemark-ET parameters with given
        argparse parser

        parser -- the argparse parser
        """
        genemark_ET_group = parser.add_argument_group('Genemark-ET parameters')
        genemark_ET_group.add_argument('--et-score', dest="et_score", type=int, help='Minimum score of intron in initiation of the ET algorithm. Default %s''' %str(self.et_score)) 
        genemark_ET_group.add_argument('--additional-genemark-ET-options', dest="additional_genemark_ET_options", help='Additional genemark-ET options to run it, see genemark documentation for more information.')
        genemark_ET_group.add_argument('--genemark-ET-weights', dest="genemark_ET_weights", nargs="+", help='Weights given to genemark-ET predictions when running EVM. Specify the weight for each EVM run separated by a space. Example  3 3 3 ') 

    def register_spaln(self, parser):
        """Register all spaln parameters with given
        argparse parser

        parser -- the argparse parser
        """
        spaln_group = parser.add_argument_group('Spaln parameters')
        spaln_group.add_argument('--ya', dest="ya", type=int, help='''Stringency of splice site. 0->3: strong->weak. Default %s''' % str(self.ya))
        spaln_group.add_argument('--M', dest="M", type=int, help='''Number of outputs per query (max=4 for genome vs cDNA|protein). Default %s''' % str(self.M))
        spaln_group.add_argument('--O', dest="O", type=int, help='''0:Gff3_gene; 1:alignment; 2:Gff3_match; 3:Bed; 4:exon-inf; 5:intron-inf; 6:cDNA; 7:translated; 8: block-only; 12: binary. Default %s''' % str(self.O))
        spaln_group.add_argument('--Q', dest="Q", type=int, help='''0:DP; 1-3: HSP-Search; 4-7; Block-Search. Default %s''' % str(self.Q))
        spaln_group.add_argument('--t', dest="t", type=int, help='''Number of threads. Default %s''' % str(self.t))
        spaln_group.add_argument('--spaln-weights', dest="spaln_weights", nargs="+", help='Weights given to spaln mappings when running EVM. Specify the weight for each EVM run separated by a space. Example  10 8 10 ') 
        spaln_group.add_argument('--additional-spaln-options', dest="additional_spaln_options", help='Additional spaln options to run it, see spaln help for more information about the possible options.')

    def register_pasa(self, parser):
        """Register all pasa parameters with given
        argparse parser

        parser -- the argparse parser
        """
        pasa_group = parser.add_argument_group('Pasa parameters')
        pasa_group.add_argument('--pasa-CPU', dest="pasa_CPU", type=int, help='''Number of pasa_CPUs to run Pasa. Default %s''' % str(self.pasa_CPU))
        pasa_group.add_argument('--pasa-step', dest="pasa_step", type=int, help='''Step from where to start running Pasa. Default %s''' % str(self.pasa_step))
        pasa_group.add_argument('--pasa-weights', dest="pasa_weights", nargs="+", help='Weights given to pasa mappings when running EVM. Specify the weight for each EVM run separated by a space. Example  8 10 8 ')
        pasa_group.add_argument('--create-database', dest="create_database", action="store_true", help='If specified, create pasa database.')

    def register_transdecoder(self, parser):
        """Register all transdecoder parameters with given
        argparse parser

        parser -- the argparse parser
        """
        transdecoder_group = parser.add_argument_group('Transdecoder parameters')
        transdecoder_group.add_argument('--transdecoder-weights', dest="transdecoder_weights", nargs="+", help='Weights given to pasa transdecodergff3 output file  when running EVM. Specify the weight for each EVM run separated by a space. Example 3 2 3')

    def register_evm(self, parser):
        """Register all evm parameters with given
        argparse parser

        parser -- the argparse parser
        """
        evm_group = parser.add_argument_group('Evm parameters')
        evm_group.add_argument('--evm-script', dest="evm_script", help='Script to run EVM. Default %s''' % str(self.evm_script))
        
    def check_parameters(self,args):
        """Check parameters consistency
            
           args -- set of parsed arguments"""

        if args.pipeline_HOME:
            self.pipeline_HOME = args.pipeline_HOME

        if args.spaln_weights:
            self.spaln_weights = args.spaln_weights

        if args.run_geneid or args.run_genemark or args.run_glimmer or args.run_augustus or args.run_geneid_introns or args.run_augustus_introns:
            if args.genome_masked == None:
                print "Sorry! No masked genome fasta file defined"
            else:
                if not os.path.exists(args.genome_masked):
                    print args.genome_masked + " not found" 
                else:
                    args.genome_masked = os.path.abspath(args.genome_masked)
        
        if args.run_pasa or args.run_transdecoder or args.run_evm or args.run_update:
            if args.genome == None:
                print "Sorry! No genome fasta file defined"
            else:
                if not os.path.exists(args.genome):
                    print args.genome + " not found" 
                else:
                    args.genome = os.path.abspath(args.genome)


        if args.annotation_step or args.annotation_version:
            if args.annotation_step:
                self.annotation_step = args.annotation_step
            if args.annotation_version:
                self.annotation_version = args.annotation_version

      
            if args.EVM_inputs:
                self.EVM_inputs = args.EVM_inputs
            else:
                self.EVM_inputs = "step0" + str(int(self.annotation_step) + 1) + "_EVM.V" + self.annotation_version 

            if args.output_dir:
                self.output_dir = args.output_dir
            else:
                self.output_dir = "step0" + str(self.annotation_step) + "_annotation_pipeline.V" + self.annotation_version 
      
            if args.update_dir:
                self.update_dir = args.update_dir
            else:
                self.update_dir = "step0" + str(int(self.annotation_step) + 2) + "_annotation_update.V" + self.annotation_version  + "/" 

            if args.augustus_prediction:
                self.augustus_prediction = args.augustus_prediction
            else: 
                self.augustus_prediction = self.output_dir + "/gene_predictions/augustus/augustus_gene_prediction.gff3"

            if args.augustus_preEVM:
                self.augustus_preEVM = args.augustus_preEVM
            else: 
                self.augustus_preEVM = self.output_dir + "/gene_predictions/augustus/augustus_preEVM.gff3"

            if args.geneid_prediction:
                self.geneid_prediction = args.geneid_prediction
            else: 
                self.geneid_prediction = self.output_dir + "/gene_predictions/geneid/geneid.gff3"

            if args.geneid_preEVM:
                self.geneid_preEVM = args.geneid_preEVM
            else: 
                self.geneid_preEVM = self.output_dir + "/gene_predictions/geneid/geneid_preEVM.gff3"

            if args.genemark_prediction:
                self.genemark_prediction = args.genemark_prediction
            else: 
                self.genemark_prediction = self.output_dir + "/gene_predictions/genemark.gtf"

            if args.genemark_preEVM:
                self.genemark_preEVM = args.genemark_preEVM
            else: 
                self.genemark_preEVM = self.output_dir + "/gene_predictions/genemark_preEVM.gff3"
 
            if args.glimmer_prediction:
                self.glimmer_prediction = args.glimmer_prediction
            else: 
                self.glimmer_prediction = self.output_dir + "/gene_predictions/glimmer.gff3" 

            if args.glimmer_preEVM:
                self.glimmer_preEVM = args.glimmer_preEVM
            else: 
                self.glimmer_preEVM = self.output_dir + "/gene_predictions/glimmer_preEVM.gff3"

            if args.geneid_introns_prediction:
                self.geneid_introns_prediction = args.geneid_introns_prediction
            else: 
                self.geneid_introns_prediction = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns.gff3"

            if args.geneid_introns_preEVM:
                self.geneid_introns_preEVM = args.geneid_introns_preEVM
            else: 
                self.geneid_introns_preEVM = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns_preEVM.gff3"

            if args.augustus_introns_prediction:
                self.augustus_introns_prediction = args.augustus_introns_prediction
            else: 
                self.augustus_introns_prediction = self.output_dir + "/gene_predictions/augustus_with_introns/augustus_introns.gff3"

            if args.augustus_introns_preEVM:
                self.augustus_introns_preEVM = args.augustus_introns_preEVM
            else: 
                self.augustus_introns_preEVM = self.output_dir + "/gene_predictions/augustus_with_introns/augustus_introns_preEVM.gff3"

            if args.genemark_ET_prediction:
                self.genemark_ET_prediction = args.genemark_ET_prediction
            else: 
                self.genemark_ET_prediction = self.output_dir + "/gene_predictions/genemark-ET.gtf"

            if args.genemark_ET_preEVM:
                self.genemark_ET_preEVM = args.genemark_ET_preEVM
            else: 
                self.genemark_ET_preEVM = self.output_dir + "/gene_predictions/genemark-ET_preEVM.gff3"

            if args.spaln_gene:
                self.spaln_gene = args.spaln_gene
            else: 
                self.spaln_gene =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_gene.gff3" 

            if args.spaln_cds:
                self.spaln_cds = args.spaln_cds
            else: 
                self.spaln_cds =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_cds.gff3" 

            if args.pasa_dir:
                self.pasa_dir = args.pasa_dir
            else: 
                self.pasa_dir = self.output_dir + "/protein_and_transcript_mappings/pasa/" 

        else:
            if args.output_dir:
                self.output_dir = args.output_dir

                if args.augustus_prediction:
                    self.augustus_prediction = args.augustus_prediction
                else: 
                    self.augustus_prediction = self.output_dir + "/gene_predictions/augustus/augustus_gene_prediction.gff3"
                
                if args.augustus_preEVM:
                    self.augustus_preEVM = args.augustus_preEVM
                else: 
                    self.augustus_preEVM = self.output_dir + "/gene_predictions/augustus/augustus_preEVM.gff3"

                if args.geneid_prediction:
                    self.geneid_prediction = args.geneid_prediction
                else: 
                    self.geneid_prediction = self.output_dir + "/gene_predictions/geneid/geneid.gff3"

                if args.geneid_preEVM:
                    self.geneid_preEVM = args.geneid_preEVM
                else: 
                    self.geneid_preEVM = self.output_dir + "/gene_predictions/geneid/geneid_preEVM.gff3"

                if args.genemark_preEVM:
                    self.genemark_preEVM = args.genemark_preEVM
                else: 
                    self.genemark_preEVM = self.output_dir + "/gene_predictions/genemark_preEVM.gff3"

                if args.glimmer_prediction:
                    self.glimmer_prediction = args.glimmer_prediction
                else: 
                    self.glimmer_prediction = self.output_dir + "/gene_predictions/glimmer.gff3" 

                if args.glimmer_preEVM:
                    self.glimmer_preEVM = args.glimmer_preEVM
                else: 
                    self.glimmer_preEVM = self.output_dir + "/gene_predictions/glimmer_preEVM.gff3"

                if args.geneid_introns_prediction:
                    self.geneid_introns_prediction = args.geneid_introns_prediction
                else: 
                    self.geneid_introns_prediction = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns.gff3"

                if args.geneid_introns_preEVM:
                    self.geneid_introns_preEVM = args.geneid_introns_preEVM
                else: 
                    self.geneid_introns_preEVM = self.output_dir + "/gene_predictions/geneid_with_introns/geneid_introns_preEVM.gff3"

                if args.augustus_introns_prediction:
                    self.augustus_introns_prediction = args.augustus_introns_prediction
                else: 
                    self.augustus_introns_prediction = self.output_dir + "/gene_predictions/augustus_with_introns/augustus_introns.gff3"

                if args.augustus_introns_preEVM:
                    self.augustus_introns_preEVM = args.augustus_introns_preEVM
                else: 
                    self.augustus_introns_preEVM = self.output_dir + "/gene_predictions/augustus_with_introns/augustus_introns_preEVM.gff3"

                if args.genemark_ET_prediction:
                    self.genemark_ET_prediction = args.genemark_ET_prediction
                else: 
                    self.genemark_ET_prediction = self.output_dir + "/gene_predictions/genemark-ET.gtf"

                if args.genemark_ET_preEVM:
                    self.genemark_ET_preEVM = args.genemark_ET_preEVM
                else: 
                    self.genemark_ET_preEVM = self.output_dir + "/gene_predictions/genemark-ET_preEVM.gff3"

                if args.spaln_gene:
                    self.spaln_gene = args.spaln_gene
                else: 
                    self.spaln_gene =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_gene.gff3" 

                if args.spaln_cds:
                    self.spaln_cds = args.spaln_cds
                else: 
                    self.spaln_cds =  self.output_dir + "/protein_and_transcript_mappings/spaln/proteins_spaln_cds.gff3" 
              
                if args.pasa_dir:
                    self.pasa_dir = args.pasa_dir
                else: 
                    self.pasa_dir = self.output_dir + "/protein_and_transcript_mappings/pasa/"             

            else:
                if args.augustus_prediction:
                    self.augustus_prediction = args.augustus_prediction
                
                if args.augustus_preEVM:
                    self.augustus_preEVM = args.augustus_preEVM

                if args.geneid_prediction:
                    self.geneid_prediction = args.geneid_prediction

                if args.geneid_preEVM:
                    self.geneid_preEVM = args.geneid_preEVM

                if args.genemark_prediction:
                    self.genemark_prediction = args.genemark_prediction

                if args.genemark_preEVM:
                    self.genemark_preEVM = args.genemark_preEVM

                if args.glimmer_prediction:
                    self.glimmer_prediction = args.glimmer_prediction

                if args.glimmer_preEVM:
                    self.glimmer_preEVM = args.glimmer_preEVM

                if args.geneid_introns_prediction:
                    self.geneid_introns_prediction = args.geneid_introns_prediction

                if args.geneid_introns_preEVM:
                    self.geneid_introns_preEVM = args.geneid_introns_preEVM

                if args.augustus_introns_prediction:
                    self.augustus_introns_prediction = args.augustus_introns_prediction

                if args.augustus_introns_preEVM:
                    self.augustus_introns_preEVM = args.augustus_introns_preEVM

                if args.genemark_ET_prediction:
                    self.genemark_ET_prediction = args.genemark_ET_prediction

                if args.genemark_ET_preEVM:
                    self.genemark_ET_preEVM = args.genemark_ET_preEVM
          
                if args.spaln_gene:
                    self.spaln_gene = args.spaln_gene

                if args.spaln_cds:
                    self.spaln_cds = args.spaln_cds
   
                if args.pasa_dir:
                    self.pasa_dir = args.pasa_dir

            if args.EVM_inputs:
                self.EVM_inputs = args.EVM_inputs

            if args.update_dir:
                self.update_dir = args.update_dir

        if args.run_augustus or args.run_augustus_introns:
            if args.species == None:
                print "Sorry! No species for augustus defined"  
             
            if args.sample:
                self.sample = args.sample

            if args.min_intron_len:
                self.min_intron_len = args.min_intron_len

            if args.alternatives_from_sampling:
                if args.alternatives_from_sampling == "false" or args.alternatives_from_sampling == "true":
                    self.alternatives_from_sampling = args.alternatives_from_sampling
                else:
                    print "Sorry! --alternatives-from-sampling must be true or false!"

            if args.uniqueGeneId:
                if args.uniqueGeneId == "false" or args.uniqueGeneId == "true":
                    self.uniqueGeneId = args.uniqueGeneId
                else:
                    print "Sorry! --uniqueGeneId must be true or false!"

            if args.gff3:
                if args.gff3 == "ON" or args.gff3 == "OFF":
                    self.gff3 = args.gff3
                else:
                    print "Sorry! --gff3 must be ON or OFF!"
 
            if args.noInFrameStop:
                if args.noInFrameStop == "false" or args.noInFrameStop == "true":
                    self.noInFrameStop = args.noInFrameStop
                else:
                    print "Sorry! --noInFrameStop must be true or false!"

            if args.singlestrand:
                if args.singlestrand == "false" or args.singlestrand == "true":
                    self.singlestrand = args.singlestrand
                else:
                    print "Sorry! --singlestrand must be true or false!"

            if args.strand:
                if args.strand == "both" or args.strand == "forward" or args.strand == "backward":
                    self.strand = args.strand
                else:
                    print "Sorry! --strand must be forward, backward or both!"

        if args.run_augustus:
            if args.augustus_weights:
                self.augustus_weights = args.augustus_weights

            if args.additional_augustus_options:
                self.additional_augustus_options = args.additional_augustus_options


        if args.run_geneid or args.run_update or args.run_geneid_introns:
            if args.geneid_parameters == None:
                print "Sorry! No geneid parameters file defined" 
            else:
                if not os.path.exists(args.geneid_parameters):
                    print args.geneid_parameters + " not found" 
                else:
                    args.geneid_parameters = os.path.abspath(args.geneid_parameters) 

        if args.run_geneid:
            if args.geneid_weights:
                self.geneid_weights = args.geneid_weights   

            if args.geneid_options:
                self.geneid_options = args.geneid_options

        if args.run_genemark or args.run_genemark_ET:
 
            if args.min_contig:
                self.min_contig = args.min_contig

            if args.max_contig:
                self.max_contig = args.max_contig

            if args.max_gap:
                self.max_gap = args.max_gap

            if args.genemark_cores:
                self.genemark_cores = args.genemark_cores

        if args.run_genemark:      
            if args.genemark_weights:
                self.genemark_weights = args.genemark_weights

            if args.additional_genemark_options:
                self.additional_genemark_options = args.additional_genemark_options

        if args.run_glimmer:
            if args.glimmer_directory == None:
                print "Sorry! No glimmer trained directory given"
            else:
                if not os.path.exists(args.glimmer_directory):
                    print args.glimmer_directory + " not found" 
                else:
                    args.glimmer_directory = os.path.abspath(args.glimmer_directory)
  
            if args.glimmer_weights:
                self.glimmer_weights = args.glimmer_weights

        if args.run_geneid_introns:
            if args.geneid_introns_weights:
                self.geneid_introns_weights = args.geneid_introns_weights   

            if args.geneid_introns_options:
                self.geneid_introns_options = args.geneid_introns_options

        if args.run_geneid_introns or args.run_augustus_introns or args.run_genemark_ET:
            if args.junctions == None and args.incoding_junctions == None:
                print "Sorry! No junctions gff file given."
            else:
                if args.junctions:
                    if not os.path.exists(args.junctions):
                        print args.junctions + " not found" 
                    else:
                        args.junctions = os.path.abspath(args.junctions)
                if args.incoding_junctions:
                    if not os.path.exists(args.incoding_junctions):
                        print args.incoding_junctions + " not found" 
                    else:
                        args.incoding_junctions = os.path.abspath(args.incoding_junctions)

        if args.run_augustus_introns:
            if args.extrinsic_file_augustus_introns == None:
                print "Sorry! No extrinsic file for augustus with introns given."
            else:
                if not os.path.exists(args.extrinsic_file_augustus_introns):
                    print args.extrinsic_file_augustus_introns + " not found" 
                else:
                    args.extrinsic_file_augustus_introns = os.path.abspath(args.extrinsic_file_augustus_introns)

            if args.augustus_introns_weights:
                self.augustus_introns_weights = args.augustus_introns_weights   

            if args.additional_augustus_introns_options:
                self.additional_augustus_introns_options = args.additional_augustus_introns_options

        if args.run_genemark_ET:
            if args.et_score:
                self.et_score = args.et_score

            if args.genemark_ET_weights:
                self.genemark_ET_weights = args.genemark_ET_weights   

            if args.additional_genemark_ET_options:
                self.additional_genemark_ET_options = args.additional_genemark_ET_options


        if args.run_spaln:
            if args.proteins == None:
                print "Sorry! No protein evidence file given."
            else:
                if not os.path.exists(args.proteins):
                    print args.proteins + " not found" 
                else:
                    args.proteins = os.path.abspath(args.proteins)
          
            if args.ya:
                self.ya = args.ya

            if args.M:
                self.M = args.M
 
            if args.O:
                self.O = args.O
        
            if args.Q:
                self.Q = args.Q

            if args.t:
                self.t = args.t

            if args.additional_spaln_options:
                self.additional_spaln_options = args.additional_spaln_options
          
        if args.run_pasa or args.run_transdecoder or args.run_update:
            if args.pasadb == None:
                print "Sorry! No pasadb name given"

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

            if args.pasa_CPU:
                self.pasa_CPU = args.pasa_CPU

            if args.pasa_step:
                self.pasa_step = args.pasa_step

            if args.pasa_weights:
                self.pasa_weights = args.pasa_weights

        if args.run_transdecoder:
            if args.transdecoder_weights:
                self.transdecoder_weights = args.transdecoder_weights

        if args.cufflinks:
            if not os.path.exists(args.cufflinks):
                print args.cufflinks + " not found" 
            else:
                args.cufflinks = os.path.abspath(args.cufflinks)   

        if args.run_evm:
            if args.evm_script:
                if not os.path.exists(args.evm_script):
                    print args.evm_script + " not found" 
                else:
                    self.evm_script = os.path.abspath(args.evm_script)

        if args.run_update:
            if args.update_config == None:
                print "Sorry! No pasa update configuration file found."
            else:
                if not os.path.exists(args.update_config):
                    print args.update_config + " not found" 
                else:
                    args.update_config = os.path.abspath(args.update_config)

            if args.project_name == None:
                print "Sorry! No project name and annotation version found."
            elif (len (args.project_name) != 2):
                print "--project-name needs 2 arguments: the name of the project and the version of the annotation"

    def storeGeneralParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.generalParameters["run_geneid"] = args.run_geneid
        self.generalParameters["run_augustus"] = args.run_augustus
        self.generalParameters["run_genemark"] = args.run_genemark
        self.generalParameters["run_glimmer"] = args.run_glimmer
        self.generalParameters["run_geneid_introns"] = args.run_geneid_introns
        self.generalParameters["run_augustus_introns"] = args.run_augustus_introns
        self.generalParameters["run_genemark_ET"] = args.run_genemark_ET
        self.generalParameters["run_spaln"] = args.run_spaln
        self.generalParameters["run_pasa"] = args.run_pasa
        self.generalParameters["run_transdecoder"] = args.run_transdecoder
        self.generalParameters["run_evm"] = args.run_evm
        self.generalParameters["run_update"] = args.run_update
        self.generalParameters["pipeline_HOME"] = self.pipeline_HOME
        self.allParameters  ["Parameters"] = self.generalParameters

    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
      #  self.inputParameters["jsonFile"] = args.jsonFile
        self.inputParameters["genome"] = args.genome
        self.inputParameters["genome_masked"] = args.genome_masked
        self.inputParameters["species"] = args.species    
        self.inputParameters["geneid_parameters"] = args.geneid_parameters
        self.inputParameters["glimmer_directory"] = args.glimmer_directory
        self.inputParameters["junctions"] = args.junctions
        self.inputParameters["incoding_junctions"] = args.incoding_junctions
        self.inputParameters["extrinsic_file_augustus_introns"] = args.extrinsic_file_augustus_introns
        self.inputParameters["proteins"] = args.proteins
        self.inputParameters["pasadb"] = args.pasadb
        self.inputParameters["transcripts"] = args.transcripts
        self.inputParameters["pasa_config"] = args.pasa_config
        self.inputParameters["cufflinks"] = args.cufflinks
        self.inputParameters["update_config"] = args.update_config
        self.inputParameters["project_name"] = args.project_name
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates output parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """

        self.outputParameters["annotation_step"] = self.annotation_step
        self.outputParameters["annotation_version"] = self.annotation_version
        self.outputParameters["output_dir"] = self.output_dir
        self.outputParameters["EVM_inputs"] = self.EVM_inputs
        self.outputParameters["update_dir"] = self.update_dir
        self.outputParameters["augustus_prediction"] = self.augustus_prediction
        self.outputParameters["augustus_preEVM"] = self.augustus_preEVM
        self.outputParameters["geneid_prediction"] = self.geneid_prediction
        self.outputParameters["geneid_preEVM"] = self.geneid_preEVM
        self.outputParameters["genemark_prediction"] = self.genemark_prediction
        self.outputParameters["genemark_preEVM"] = self.genemark_preEVM
        self.outputParameters["glimmer_prediction"] = self.glimmer_prediction
        self.outputParameters["glimmer_preEVM"] = self.glimmer_preEVM
        self.outputParameters["geneid_introns_prediction"] = self.geneid_introns_prediction
        self.outputParameters["geneid_introns_preEVM"] = self.geneid_introns_preEVM
        self.outputParameters["augustus_introns_prediction"] = self.augustus_introns_prediction
        self.outputParameters["augustus_introns_preEVM"] = self.augustus_introns_preEVM
        self.outputParameters["genemark_ET_prediction"] = self.genemark_ET_prediction
        self.outputParameters["genemark_ET_preEVM"] = self.genemark_ET_preEVM
        self.outputParameters["spaln_gene"] = self.spaln_gene
        self.outputParameters["spaln_cds"] = self.spaln_cds
        self.outputParameters["pasa_dir"] = self.pasa_dir
        self.allParameters  ["Outputs"] = self.outputParameters

    def storeAugustusParameters(self,args):
        """Updates augustus parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusParameters["alternatives_from_sampling"] = self.alternatives_from_sampling
        self.augustusParameters["uniqueGeneId"] = self.uniqueGeneId
        self.augustusParameters["gff3"] = self.gff3
        self.augustusParameters["sample"] = self.sample
        self.augustusParameters["noInFrameStop"] = self.noInFrameStop
        self.augustusParameters["maxtracks"] = self.maxtracks
        self.augustusParameters["singlestrand"] = self.singlestrand
        self.augustusParameters["strand"] = self.strand
        self.augustusParameters["min_intron_len"] = self.min_intron_len
        self.augustusParameters["augustus_weights"] = self.augustus_weights
        self.augustusParameters["additional_augustus_options"] = self.additional_augustus_options      
        self.allParameters ["augustus"] = self.augustusParameters

    def storeGeneidParameters(self,args):
        """Updates geneid parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.geneidParameters["geneid_weights"] = self.geneid_weights
        self.geneidParameters["geneid_options"] = self.geneid_options       
        self.allParameters ["geneid"] = self.geneidParameters
        
    def storeGenemarkParameters(self,args):
        """Updates genemark parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genemarkParameters["max_contig"] = self.max_contig
        self.genemarkParameters["min_contig"] = self.min_contig
        self.genemarkParameters["max_gap"] = self.max_gap
        self.genemarkParameters["genemark_cores"] = self.genemark_cores
        self.genemarkParameters["genemark_weights"] = self.genemark_weights
        self.genemarkParameters["additional_genemark_options"] = self.additional_genemark_options   
        self.allParameters ["genemark"] = self.genemarkParameters

    def storeGlimmerParameters(self,args):
        """Updates glimmer parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.glimmerParameters["glimmer_weights"] = self.glimmer_weights
        self.allParameters ["glimmer"] = self.glimmerParameters

    def storeGeneidIntronsParameters(self,args):
        """Updates geneid with introns parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.geneidIntronsParameters["geneid_introns_weights"] = self.geneid_introns_weights
        self.geneidIntronsParameters["geneid_introns_options"] = self.geneid_introns_options       
        self.allParameters ["geneid_introns"] = self.geneidIntronsParameters
        
    def storeAugustusIntronsParameters(self,args):
        """Updates augustus with introns parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.augustusIntronsParameters["augustus_introns_weights"] = self.augustus_introns_weights
        self.augustusIntronsParameters["additional_augustus_introns_options"] = self.additional_augustus_introns_options       
        self.allParameters ["augustus_introns"] = self.augustusIntronsParameters

    def storeGenemarkETParameters(self,args):
        """Updates genemark parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genemarkETParameters["genemark_ET_weights"] = self.genemark_ET_weights
        self.genemarkETParameters["et_score"] = self.et_score
        self.genemarkETParameters["additional_genemark_ET_options"] = self.additional_genemark_ET_options 
        self.allParameters ["genemark-ET"] = self.genemarkETParameters

    def storeSpalnParameters(self,args):
        """Updates spaln parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.spalnParameters["ya"] = self.ya
        self.spalnParameters["ya"] = self.ya
        self.spalnParameters["M"] = self.M
        self.spalnParameters["O"] = self.O
        self.spalnParameters["Q"] = self.Q
        self.spalnParameters["t"] = self.t
        self.spalnParameters["spaln_weights"] = self.spaln_weights
        self.spalnParameters["additional_spaln_options"] = self.additional_spaln_options      
        self.allParameters ["spaln"] = self.spalnParameters

    def storePasaParameters(self,args):
        """Updates pasa parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pasaParameters["pasa_CPU"] = self.pasa_CPU
        self.pasaParameters["create_database"] = args.create_database
        self.pasaParameters["pasa_step"] = self.pasa_step
        self.pasaParameters["pasa_weights"] = self.pasa_weights
        self.allParameters["pasa"] = self.pasaParameters

    def storeTransdecoderParameters(self,args):
        """Updates transdecoder parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.transdecoderParameters["transdecoder_weights"] = self.transdecoder_weights
        self.allParameters["transdecoder"] = self.transdecoderParameters

    def storeEvmParameters(self,args):
        """Updates evm parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.evmParameters["evm_script"] = self.evm_script
        self.allParameters["evm"] = self.evmParameters

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
configManager.storeAugustusParameters(args)
configManager.storeGeneidParameters(args)
configManager.storeGenemarkParameters(args)
configManager.storeGlimmerParameters(args)
configManager.storeGeneidIntronsParameters(args)
configManager.storeAugustusIntronsParameters(args)
configManager.storeGenemarkETParameters(args)
configManager.storeSpalnParameters(args)
configManager.storePasaParameters(args)
configManager.storeTransdecoderParameters(args)
configManager.storeEvmParameters(args)

#4. Store JSON file
with open(args.jsonFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)




#Class Create Configuration File
   # Define default values

   # Parser arguments

   # store arguments to map

#Create object class Configuration File

#Call argument parsing

#Call store arguments method by pipeline step

#Create super map fro every map pipeline step

#store super to Json file
