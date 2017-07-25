#!/usr/bin/env python
import fileinput
import sys

__author__ = 'adamkoziol'

"""
Uses logic from https://github.com/katholt/srst2/blob/master/scripts/srst2.py
"""


def editheaders():
    """Edits the headers of SAM files to remove 'secondary alignments'"""
    # Read stdin - this will be the output from samtools view
    for line in fileinput.input():
        try:
            # Get the flag value from the input
            columns = line.split('\t')
            # The FLAG is in the second column
            flag = int(columns[1])
            # Subtracts 256 from the flag if the & bitwise operator evaluates to true
            # See http://www.tutorialspoint.com/python/bitwise_operators_example.htm
            # For the test case, flags of 256 became 0, and flags of 272 became 16
            columns[1] = str((flag - 256) if (flag & 256) else flag)
            # update = [columns[0], str(flag), columns[2:]]
            sys.stdout.write('\t'.join(columns))
        # Don't fail on IOErrors, or ValueErrors, and still print line to stdout
        except (IOError, ValueError):
            sys.stdout.write(line)
            pass
    # Try except statements to get rid of file closing errors
    try:
        sys.stdout.flush()
        sys.stdout.close()
    except:
        pass
    try:
        sys.stderr.close()
    except:
        pass

if __name__ == '__main__':
    # Run the script
    editheaders()
