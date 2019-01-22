#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import GenObject
from datetime import datetime
import logging
import os
__author__ = 'adamkoziol'


class Reporter(object):

    def reporter(self):
        """
        Creates the metadata report by pulling specific attributes from the metadata objects
        """
        logging.info('Creating summary report')
        header = '{}\n'.format(','.join(self.headers))
        # Create a string to store all the results
        data = str()
        for sample in self.metadata:
            # Add the value of the appropriate attribute to the results string
            data += GenObject.returnattr(sample, 'name')
            # SampleName
            data += GenObject.returnattr(sample.run, 'SamplePlate')
            # Genus
            data += GenObject.returnattr(sample.general, 'closestrefseqgenus')
            # SequencingDate
            data += GenObject.returnattr(sample.run, 'Date')
            # Analyst
            data += GenObject.returnattr(sample.run, 'InvestigatorName')
            # SamplePurity
            data += GenObject.returnattr(sample.confindr, 'contam_status')
            # N50
            n50 = GenObject.returnattr(sample.quality_features_polished, 'n50',
                                       number=True)
            if n50 != '-,':
                data += n50
            else:
                data += '0,'
            # NumContigs
            data += GenObject.returnattr(sample.quality_features_polished, 'num_contigs',
                                         number=True)
            # TotalLength
            data += GenObject.returnattr(sample.quality_features_polished, 'genome_length',
                                         number=True)
            # MeanInsertSize
            data += GenObject.returnattr(sample.mapping, 'MeanInsertSize',
                                         number=True)
            # InsertSizeSTD
            data += GenObject.returnattr(sample.mapping, 'StdInsertSize',
                                         number=True)
            # AverageCoverageDepth
            data += GenObject.returnattr(sample.mapping, 'MeanCoveragedata',
                                         number=True)
            # CoverageDepthSTD
            data += GenObject.returnattr(sample.mapping, 'StdCoveragedata',
                                         number=True)
            # PercentGC
            data += GenObject.returnattr(sample.quality_features_polished, 'gc',
                                         number=True)
            # MASH_ReferenceGenome
            data += GenObject.returnattr(sample.mash, 'closestrefseq')
            # MASH_NumMatchingHashes
            data += GenObject.returnattr(sample.mash, 'nummatches')
            # 16S_result
            data += GenObject.returnattr(sample.sixteens_full, 'sixteens_match')
            # rMLST_Result
            try:
                # If the number of matches to the closest reference profile is 53, return the profile number
                if sample.rmlst.matches == 53:
                    rmlst_seq_type = GenObject.returnattr(sample.rmlst, 'sequencetype')
                    rmlst_seq_type = rmlst_seq_type if rmlst_seq_type != 'ND,' else 'new,'
                    data += rmlst_seq_type
                else:
                    # Otherwise the profile is set to new
                    data += 'new,'
            except AttributeError:
                data += 'new,'
            # MLST_Result
            try:
                if sample.mlst.matches == 7:
                    data += GenObject.returnattr(sample.mlst, 'sequencetype')
                else:
                    data += 'new,'
                    # # Create a set of all the genes present in the results (gene name split from allele)
                    # mlst_gene_set = {gene.split('_')[0] for gene in sample.mlst.results}
                    # # If there are all the genes present, but no perfect match to a reference profile, state that
                    # # the profile is new
                    # if len(mlst_gene_set) == 7:
                    #     data += 'new,'
                    # # Otherwise indicate that the profile is ND
                    # else:
                    #     data += 'ND,'
            except AttributeError:
                data += 'new,'
            # MLST_gene_X_alleles
            try:
                # Create a set of all the genes present in the results (gene name split from allele)
                gene_set = {gene.split('_')[0] for gene in sample.mlst.results}
                for gene in sorted(gene_set):
                    allele_list = list()
                    # Determine all the alleles that are present for each gene
                    for allele in sample.mlst.results:
                        if gene in allele:
                            allele_list.append(allele)
                    # If there is more than one allele in the sample, add both to the string separated by a ';'
                    if len(allele_list) > 1:
                        data += '{},'.format(';'.join(allele_list))
                    # Otherwise add the only allele
                    else:
                        data += allele_list[0] + ','
                # If there are fewer than seven matching alleles, add a ND for each missing result
                if len(gene_set) < 7:
                    data += (7 - len(gene_set)) * 'ND,'
            except AttributeError:
                # data += '-,-,-,-,-,-,-,'
                data += 'ND,ND,ND,ND,ND,ND,ND,'
            # CoreGenesPresent
            data += GenObject.returnattr(sample.coregenome, 'coreresults')
            # E_coli_Serotype
            try:
                # If no O-type was found, set the output to be O-untypeable
                if ';'.join(sample.serosippr.o_set) == '-':
                    otype = 'O-untypeable'
                else:
                    otype = '{oset} ({opid})'.format(oset=';'.join(sample.serosippr.o_set),
                                                     opid=sample.serosippr.best_o_pid)
                # Same as above for the H-type
                if ';'.join(sample.serosippr.h_set) == '-':
                    htype = 'H-untypeable'

                else:
                    htype = '{hset} ({hpid})'.format(hset=';'.join(sample.serosippr.h_set),
                                                     hpid=sample.serosippr.best_h_pid)
                serotype = '{otype}:{htype},'.format(otype=otype,
                                                     htype=htype)
                # Add the serotype to the data string unless neither O-type not H-type were found; add ND instead
                data += serotype if serotype != 'O-untypeable:H-untypeable,' else 'ND,'
            except AttributeError:
                data += 'ND,'
            # SISTR_serovar_antigen
            data += GenObject.returnattr(sample.sistr, 'serovar_antigen').rstrip(';')
            # SISTR_serovar_cgMLST
            data += GenObject.returnattr(sample.sistr, 'serovar_cgmlst')
            # SISTR_serogroup
            data += GenObject.returnattr(sample.sistr, 'serogroup')
            # SISTR_h1
            data += GenObject.returnattr(sample.sistr, 'h1').rstrip(';')
            # SISTR_h2
            data += GenObject.returnattr(sample.sistr, 'h2').rstrip(';')
            # SISTR_serovar
            data += GenObject.returnattr(sample.sistr, 'serovar')
            # GeneSeekr_Profile
            try:
                if sample.genesippr.report_output:
                    data += ';'.join(sample.genesippr.report_output) + ','
                else:
                    data += 'ND,'
            except AttributeError:
                data += 'ND,'
            # Vtyper_Profile
            data += GenObject.returnattr(sample.legacy_vtyper, 'toxinprofile')
            # AMR_Profile and resistant/sensitive status
            if sample.resfinder_assembled.pipelineresults:
                # Profile
                for resistance, resistance_set in sorted(sample.resfinder_assembled.pipelineresults.items()):
                    data += '{res}({r_set});'.format(res=resistance.replace(',', ';'),
                                                     r_set=';'.join(sorted(list(resistance_set))))
                data += ','
                # Resistant/Sensitive
                data += 'Resistant,'
            else:
                # Profile
                data += 'ND,'
                # Resistant/Sensitive
                data += 'Sensitive,'
            # Plasmid Result'
            if sample.mobrecon.pipelineresults:
                for plasmid, details in sorted(sample.mobrecon.pipelineresults.items()):
                    data += '{plasmid}({details});'.format(plasmid=plasmid,
                                                           details=details)
                data += ','
            else:
                data += 'ND,'
            # TotalPredictedGenes
            data += GenObject.returnattr(sample.prodigal, 'predictedgenestotal',
                                         number=True)
            # PredictedGenesOver3000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover3000bp',
                                         number=True)
            # PredictedGenesOver1000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover1000bp',
                                         number=True)
            # PredictedGenesOver500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover500bp',
                                         number=True)
            # PredictedGenesUnder500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesunder500bp',
                                         number=True)
            # NumClustersPF
            data += GenObject.returnattr(sample.run, 'NumberofClustersPF')
            # Percent of reads mapping to PhiX control
            data += GenObject.returnattr(sample.run, 'phix_aligned')
            # Error rate calculated from PhiX control
            data += GenObject.returnattr(sample.run, 'error_rate')
            # LengthForwardRead
            data += GenObject.returnattr(sample.run, 'forwardlength',
                                         number=True)
            # LengthReverseRead
            data += GenObject.returnattr(sample.run, 'reverselength',
                                         number=True)
            # Real time strain
            data += GenObject.returnattr(sample.run, 'Description')
            # Flowcell
            data += GenObject.returnattr(sample.run, 'flowcell')
            # MachineName
            data += GenObject.returnattr(sample.run, 'instrument')
            # PipelineVersion
            data += self.commit + ','
            # AssemblyDate
            data += datetime.now().strftime('%Y-%m-%d')
            # Append a new line to the end of the results for this sample
            data += '\n'
        # Replace any NA values with -
        cleandata = data.replace('NA', 'ND')
        with open(os.path.join(self.reportpath, 'combinedMetadata.csv'), 'w') as metadatareport:
            metadatareport.write(header)
            metadatareport.write(cleandata)

    def legacy_reporter(self):
        """
        Creates an output that is compatible with the legacy metadata reports. This method will be removed once
        a new database scheme is implemented
        """
        from collections import OrderedDict
        logging.info('Creating legacy summary report')
        row = ''
        # Create a dictionary of tuples to be printed in the final report
        for sample in self.metadata:
            data = OrderedDict([
                ('SampleName', sample.name),
                ('N50', str(sample.quality_features_polished.n50)),
                ('NumContigs', str(sample.quality_features_polished.num_contigs)),
                ('TotalLength', str(sample.quality_features_polished.genome_length)),
                ('MeanInsertSize', sample.mapping.MeanInsertSize),
                ('AverageCoverageDepth', sample.mapping.MeanCoveragedata.split("X")[0]),
                ('ReferenceGenome', sample.mash.closestrefseq),
                ('RefGenomeAlleleMatches', '-'),
                ('16sPhylogeny', sample.sixteens_full.genus),
                ('rMLSTsequenceType', sample.rmlst.sequencetype),
                ('MLSTsequencetype', sample.mlst.sequencetype),
                ('MLSTmatches', str(sample.mlst.matchestosequencetype)),
                ('coreGenome', GenObject.returnattr(sample.coregenome, 'coreresults').rstrip(',')),
                ('SeroType', '{oset}:{hset}'
                    .format(oset=';'.join(sample.serosippr.o_set),
                            hset=';'.join(sample.serosippr.h_set))),
                ('geneSeekrProfile', ';'.join(result for result, pid in sorted(sample.genesippr.results.items()))),
                ('vtyperProfile', ';'.join(sorted(sample.legacy_vtyper.toxinprofile))),
                ('percentGC', str(sample.quality_features_polished.gc)),
                ('TotalPredictedGenes', str(sample.prodigal.predictedgenestotal)),
                ('predictedgenesover3000bp', str(sample.prodigal.predictedgenesover3000bp)),
                ('predictedgenesover1000bp', str(sample.prodigal.predictedgenesover1000bp)),
                ('predictedgenesover500bp', str(sample.prodigal.predictedgenesover500bp)),
                ('predictedgenesunder500bp', str(sample.prodigal.predictedgenesunder500bp)),
                ('SequencingDate', sample.run.Date),
                ('Investigator', sample.run.InvestigatorName),
                ('TotalClustersinRun', str(sample.run.TotalClustersinRun)),
                ('NumberofClustersPF', str(sample.run.NumberofClustersPF)),
                ('PercentOfClusters', str(sample.run.PercentOfClusters)),
                ('LengthofForwardRead', str(sample.run.forwardlength)),
                ('LengthofReverseRead', str(sample.run.reverselength)),
                ('Project', str(sample.run.SampleProject)),
                ('PipelineVersion', self.commit)
            ])

            if not row:
                row += ','.join([key for key, value in data.items()])
            row += '\n'
            row += ','.join([value for key, value in data.items()])
        cleanrow = row.replace('NA', '').replace(',-,', ',,')
        with open(os.path.join(self.reportpath, 'legacy_combinedMetadata.csv'), 'w') as metadatareport:
            metadatareport.write(cleanrow)

    def clean_object(self):
        for sample in self.metadata:
            try:
                delattr(sample.coregenome, 'targetnames')
            except AttributeError:
                pass
            try:
                delattr(sample.coregenome, 'targets')
            except AttributeError:
                pass

    def __init__(self, inputobject, legacy=False):
        self.metadata = inputobject.runmetadata.samples
        self.commit = inputobject.commit
        self.reportpath = inputobject.reportpath
        self.starttime = inputobject.starttime
        # Define the headers to be used in the metadata report
        # 'AssemblyQuality',
        self.headers = ['SeqID', 'SampleName', 'Genus', 'SequencingDate', 'Analyst', 'SamplePurity',
                        'N50', 'NumContigs', 'TotalLength', 'MeanInsertSize', 'InsertSizeSTD',
                        'AverageCoverageDepth', 'CoverageDepthSTD', 'PercentGC', 'MASH_ReferenceGenome',
                        'MASH_NumMatchingHashes', '16S_result', 'rMLST_Result', 'MLST_Result', 'MLST_gene_1_allele',
                        'MLST_gene_2_allele', 'MLST_gene_3_allele', 'MLST_gene_4_allele', 'MLST_gene_5_allele',
                        'MLST_gene_6_allele', 'MLST_gene_7_allele', 'CoreGenesPresent', 'E_coli_Serotype',
                        'SISTR_serovar_antigen', 'SISTR_serovar_cgMLST', 'SISTR_serogroup', 'SISTR_h1', 'SISTR_h2',
                        'SISTR_serovar', 'GeneSeekr_Profile', 'Vtyper_Profile', 'AMR_Profile',
                        'AMR Resistant/Sensitive', 'PlasmidProfile', 'TotalPredictedGenes', 'PredictedGenesOver3000bp',
                        'PredictedGenesOver1000bp', 'PredictedGenesOver500bp', "PredictedGenesUnder500bp",
                        'NumClustersPF', 'PercentReadsPhiX', 'ErrorRate', 'LengthForwardRead', 'LengthReverseRead',
                        'RealTimeStrain', 'Flowcell', 'MachineName', 'PipelineVersion', 'AssemblyDate']
        self.reporter()
        if legacy:
            self.legacy_reporter()
        # Create a database to store all the metadata
        self.clean_object()
