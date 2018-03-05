#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import GenObject, printtime
from datetime import datetime
import os
__author__ = 'adamkoziol'


class Reporter(object):

    def reporter(self):
        """
        Creates the metadata report by pulling specific attributes from the metadata objects
        """
        printtime('Creating summary report', self.starttime)
        header = '{}\n'.format(','.join(self.headers))
        # Create a string to store all the results
        data = str()
        for sample in self.metadata:
            # Add the value of the appropriate attribute to the results string
            data += GenObject.returnattr(sample, 'name')
            # SampleName
            data += GenObject.returnattr(sample.run, 'SampleName')
            # Genus
            data += GenObject.returnattr(sample.sixteens_full, 'genus')
            # SequencingDate
            data += GenObject.returnattr(sample.run, 'Date')
            # Analyst
            data += GenObject.returnattr(sample.run, 'InvestigatorName')
            # SamplePurity
            data += GenObject.returnattr(sample.confinder, 'contam_status')
            # GenomeQAML prediction
            data += GenObject.returnattr(sample.GenomeQAML, 'prediction')
            # N50
            data += GenObject.returnattr(sample.quality_features, 'n50')
            # NumContigs
            data += GenObject.returnattr(sample.quality_features, 'num_contigs')
            # TotalLength
            data += GenObject.returnattr(sample.quality_features, 'genome_length')
            # MeanInsertSize
            data += GenObject.returnattr(sample.mapping, 'MeanInsertSize')
            # InsertSizeSTD
            data += GenObject.returnattr(sample.mapping, 'StdInsertSize')
            # AverageCoverageDepth
            data += GenObject.returnattr(sample.mapping, 'MeanCoveragedata')
            # CoverageDepthSTD
            data += GenObject.returnattr(sample.mapping, 'StdCoveragedata')
            # PercentGC
            data += GenObject.returnattr(sample.quality_features, 'gc')
            # MASH_ReferenceGenome
            data += GenObject.returnattr(sample.mash, 'closestrefseq')
            # MASH_NumMatchingHashes
            data += GenObject.returnattr(sample.mash, 'nummatches')
            # 16S_result
            data += GenObject.returnattr(sample.sixteens_full, 'sixteens_match')
            # rMLST_Result
            try:
                data += GenObject.returnattr(sample.rmlst, 'sequencetype') if sample.rmlst.matches == 53 else '-,'
            except KeyError:
                data += '-,'
            # MLST_Result
            try:
                data += GenObject.returnattr(sample.mlst, 'sequencetype') if sample.mlst.matches == 7 else '-,'
            except KeyError:
                data += '-,'
            # MLST_gene_X_alleles
            try:
                for allele in sorted(sample.mlst.results):
                    data += allele + ','
                # If there are fewer than seven matching alleles, add a - for each missing result
                if len(sample.mlst.results) < 7:
                    data += (7 - len(sample.mlst.results)) * '-,'
            except KeyError:
                data += '-,-,-,-,-,-,-,'
            # CoreGenesPresent
            data += GenObject.returnattr(sample.coregenome, 'coreresults')
            # E_coli_Serotype
            try:
                serotype = '{oset} ({opid}):{hset} ({hpid}),'\
                    .format(oset=';'.join(sample.serosippr.o_set),
                            opid=sample.serosippr.best_o_pid,
                            hset=';'.join(sample.serosippr.h_set),
                            hpid=sample.serosippr.best_h_pid)
                # Make sure that the string was populated with values rather than 'NA' or '-'
                if serotype == '- (-):- (-),':
                    data += '-,'
                else:
                    data += serotype
            except KeyError:
                data += '-,'
            # SISTR_serovar_antigen
            data += GenObject.returnattr(sample.sistr, 'serovar_antigen')
            # SISTR_serovar_cgMLST
            data += GenObject.returnattr(sample.sistr, 'serovar_cgmlst')
            # SISTR_serogroup
            data += GenObject.returnattr(sample.sistr, 'serogroup')
            # SISTR_h1
            data += GenObject.returnattr(sample.sistr, 'h1') .replace(',', ';') + ','
            # SISTR_h2
            data += GenObject.returnattr(sample.sistr, 'h2').replace(',', ';') + ','
            # SISTR_serovar
            data += GenObject.returnattr(sample.sistr, 'serovar')
            # GeneSeekr_Profile
            try:
                results = list()
                for result, percentid in sorted(sample.genesippr.results.items()):
                    if float(percentid) > 70:
                        results.append(result)
                data += ';'.join(results) + ','
            except KeyError:
                data += '-,'
            # Vtyper_Profile
            try:
                data += ';'.join(sorted(sample.vtyper.profile)) + ','
            except KeyError:
                data += '-,'
            # AMR_Profile
            try:
                data += ';'.join(sorted(sample.resfinder.pipelineresults)) + ','
            except KeyError:
                data += '-,'
            # Plasmid Result'
            try:
                data += ';'.join(sorted(sample.plasmidextractor.plasmids)) + ','
            except KeyError:
                data += '-,'
            # TotalPredictedGenes
            data += GenObject.returnattr(sample.prodigal, 'predictedgenestotal')
            # PredictedGenesOver3000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover3000bp')
            # PredictedGenesOver1000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover1000bp')
            # PredictedGenesOver500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover500bp')
            # PredictedGenesUnder500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesunder500bp')
            # NumClustersPF
            data += GenObject.returnattr(sample.run, 'NumberofClustersPF')
            # LengthForwardRead
            data += GenObject.returnattr(sample.run, 'forwardlength')
            # LengthReverseRead
            data += GenObject.returnattr(sample.run, 'reverselength')
            # OtherNames
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
        cleandata = data.replace('NA', '-')
        with open(os.path.join(self.reportpath, 'combinedMetadata.csv'), 'w') as metadatareport:
            metadatareport.write(header)
            metadatareport.write(cleandata)

    def legacy_reporter(self):
        """
        Creates an output that is compatible with the legacy metadata reports. This method will be removed once
        a new database scheme is implemented
        """
        from collections import OrderedDict
        printtime('Creating legacy summary report', self.starttime)
        row = ''
        # Create a dictionary of tuples to be printed in the final report
        for sample in self.metadata:
            data = OrderedDict([
                ('SampleName', sample.name),
                ('N50', str(sample.quality_features.n50)),
                ('NumContigs', str(sample.quality_features.num_contigs)),
                ('TotalLength', str(sample.quality_features.genome_length)),
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
                ('vtyperProfile', ';'.join(sorted(sample.vtyper.profile))),
                ('percentGC', str(sample.quality_features.gc)),
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
        cleanrow = row.replace('NA', '-')
        with open(os.path.join(self.reportpath, 'legacy_combinedMetadata.csv'), 'w') as metadatareport:
            metadatareport.write(cleanrow)

    def database(self):
        """
        Enters all the metadata into a database
        """
        import sqlite3
        try:
            os.remove('{}/metadatabase.sqlite'.format(self.reportpath))
        except OSError:
            pass
        # Set the name of the database
        db = sqlite3.connect('{}/metadatabase.sqlite'.format(self.reportpath))
        # Create a cursor to allow access to the database
        cursor = db.cursor()
        # Set up the db
        cursor.execute('''
          CREATE TABLE IF NOT EXISTS Samples (
            id     INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name   TEXT UNIQUE
          )

        ''')
        # Create a variable to store the names of the header values for each individual table
        # This will store a set of all the headers from all the strains, as there can be some variability present, as
        # not all analyses are available for all taxonomic groups
        columns = dict()
        for sample in self.metadata:
            # Insert each strain name into the Samples table
            cursor.execute('''
              INSERT OR IGNORE INTO Samples (name)
              VALUES ( ? )
            ''', (sample.name, ))
            # Each header in the .json file represents a major category e.g. ARMI, geneseekr, commands, etc. and
            # will be made into a separate table
            for header in sample.datastore.items():
                # Set the table name
                tablename = header[0].replace('.', '_')
                # Allow for certain analyses, such as core genome, not being performed on all strains
                try:
                    # Create the table (if it doesn't already exist)
                    # sample_id INTEGER
                    cursor.execute('''
                      CREATE TABLE IF NOT EXISTS {} (
                      sample_id INTEGER,
                      FOREIGN KEY(sample_id) REFERENCES Samples(id)
                      )
                      '''.format(tablename))
                    # Key and value: data description and data value e.g. targets present: 1012, etc.
                    for key, value in header[1].datastore.items():
                        # Add the data header to the dictionary
                        # Clean the column names so there are no issues entering names into the database
                        cleanedcolumn = self.columnclean(key)
                        try:
                            columns[tablename].add(str(cleanedcolumn))
                        # Initialise the dictionary the first time a table name is encountered
                        except KeyError:
                            columns[tablename] = set()
                            columns[tablename].add(str(cleanedcolumn))
                except (AttributeError, IndexError):
                    pass
        # Iterate through the dictionary containing all the data headers
        for table, setofheaders in columns.items():
            # Each header will be used as a column in the appropriate table
            for cleanedcolumn in setofheaders:
                # Alter the table by adding each header as a column
                cursor.execute('''
                  ALTER TABLE {}
                  ADD COLUMN {} TEXT
                '''.format(table, cleanedcolumn))
            # Iterate through the samples and pull out the data for each table/column
            for sample in self.metadata:
                # Find the id associated with each sample in the Sample table
                cursor.execute('''
                  SELECT id from Samples WHERE name=?
                ''', (sample.name,))
                sampleid = cursor.fetchone()[0]
                # Add the sample_id to the table
                cursor.execute('''
                  INSERT OR IGNORE INTO {}
                  (sample_id) VALUES ("{}")
                 '''.format(table, sampleid))
                # Add the data to the table
                try:
                    # Find the data for each table/column
                    for item in sample[table].datastore.items():
                        # Clean the names
                        cleanedcolumn = self.columnclean(str(item[0]))
                        # Add the data to the column of the appropriate table,
                        # where the sample_id matches the current strain
                        cursor.execute('''
                          UPDATE {}
                          SET {} = ?
                          WHERE sample_id = {}
                          '''.format(table, cleanedcolumn, sampleid), (str(item[1]), ))
                except KeyError:
                    pass
        # Commit the changes to the database
        db.commit()

    @staticmethod
    def columnclean(column):
        """
        Modifies column header format to be importable into a database
        :param column: raw column header
        :return: cleanedcolumn: reformatted column header
        """
        cleanedcolumn = str(column) \
            .replace('%', 'percent') \
            .replace('(', '_') \
            .replace(')', '') \
            .replace('As', 'Adenosines') \
            .replace('Cs', 'Cytosines') \
            .replace('Gs', 'Guanines') \
            .replace('Ts', 'Thymines') \
            .replace('Ns', 'Unknowns') \
            .replace('index', 'adapterIndex')
        return cleanedcolumn

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.commit = inputobject.commit
        self.reportpath = inputobject.reportpath
        self.starttime = inputobject.starttime
        # Define the headers to be used in the metadata report
        self.headers = ['SeqID', 'SampleName', 'Genus', 'SequencingDate', 'Analyst', 'SamplePurity',
                        'AssemblyQuality', 'N50', 'NumContigs', 'TotalLength', 'MeanInsertSize', 'InsertSizeSTD',
                        'AverageCoverageDepth', 'CoverageDepthSTD', 'PercentGC', 'MASH_ReferenceGenome',
                        'MASH_NumMatchingHashes', '16S_result', 'rMLST_Result', 'MLST_Result', 'MLST_gene_1_allele',
                        'MLST_gene_2_allele', 'MLST_gene_3_allele', 'MLST_gene_4_allele', 'MLST_gene_5_allele',
                        'MLST_gene_6_allele', 'MLST_gene_7_allele', 'CoreGenesPresent', 'E_coli_Serotype',
                        'SISTR_serovar_antigen', 'SISTR_serovar_cgMLST', 'SISTR_serogroup', 'SISTR_h1', 'SISTR_h2',
                        'SISTR_serovar', 'GeneSeekr_Profile', 'Vtyper_Profile', 'AMR_Profile', 'PlasmidProfile',
                        'TotalPredictedGenes', 'PredictedGenesOver3000bp', 'PredictedGenesOver1000bp',
                        'PredictedGenesOver500bp', "PredictedGenesUnder500bp", 'NumClustersPF', 'LengthForwardRead',
                        'LengthReverseRead', 'OtherNames', 'Flowcell', 'MachineName', 'PipelineVersion', 'AssemblyDate']
        self.reporter()
        self.legacy_reporter()
        # Create a database to store all the metadata
        self.database()
