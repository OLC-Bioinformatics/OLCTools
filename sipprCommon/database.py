#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import *
__author__ = 'adamkoziol'


class Database(object):

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
            # Create a metadata object to store the new tables
            data = MetadataObject()
            data.name = sample.name
            # Insert each strain name into the Samples table
            cursor.execute('''
              INSERT OR IGNORE INTO Samples (name)
              VALUES ( ? )
            ''', (sample.name, ))
            # Each header in the .json file represents a major category e.g. ARMI, GeneSeekr, commands, etc. and
            # will be made into a separate table
            for header in sample.datastore.items():
                # Allow for certain analyses, such as core genome, not being performed on all strains
                try:

                    # Key and value: data description and data value e.g. targets present: 1012, etc.
                    for key, value in sorted(header[1].datastore.items()):
                        # Only the values consisting of dictionaries are of interest
                        if type(value) == dict:
                            # Clean the column names so there are no issues entering names into the database
                            cleanedcolumn = self.columnclean(key)
                            # Set the table name
                            tablename = '{}_{}'.format(header[0].replace('.', '_'), cleanedcolumn)
                            # Create the table (if it doesn't already exist)
                            cursor.execute('''
                                          CREATE TABLE IF NOT EXISTS {} (
                                            sample_id INTEGER
                                          )
                                          '''.format(tablename))
                            # Add the attributes with the dictionaries (values) to the metadata object
                            setattr(data, tablename, GenObject(value))
                            for gene, result in sorted(value.items()):
                                # Add the data header to the dictionary
                                try:
                                    columns[tablename].add(gene)
                                # Initialise the dictionary the first time a table name is encountered
                                except KeyError:
                                    columns[tablename] = set()
                                    columns[tablename].add(str(gene))
                except (AttributeError, IndexError):
                    pass
            self.tabledata.append(data)
        # Iterate through the dictionary containing all the data headers
        for table, setofheaders in sorted(columns.items()):
            # Each header will be used as a column in the appropriate table
            for cleanedcolumn in sorted(setofheaders):
                # Alter the table by adding each header as a column
                cursor.execute('''
                  ALTER TABLE {}
                  ADD COLUMN {} TEXT
                '''.format(table, cleanedcolumn))
            # Iterate through the samples and pull out the data for each table/column
            # for sample in self.metadata:
            for sample in self.tabledata:
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
                    for item in sorted(sample[table].datastore.items()):
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
        self.tabledata = list()
        # Create a database to store all the metadata
        self.database()
