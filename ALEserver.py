import sys, os
import subprocess
import commands
import time


__version__ = '0.01'
__versionDate__ = '15 Feb 2011'
__author__ = 'Scott Clark <sc932 at cornell dot edu>'
__description__ = 'an alignment likelihood engine'
__usage__ = """python ALEserver.py [options] <assembly file> <read file/directory> <output directory>

where basic options are:
  -h : show brief help on version and usage
"""
__fullUsage__ = """
# ALE - """ + __description__ + """
# ALEserver """ + __version__ + ' (' + __versionDate__ + """)
# Copyright (C) 2010,2011 """ + __author__ + """

""" + __usage__ + """
Input options (read left to right):
  -readDir   : Instead of a read file, a directory of read files is given
  -split     : The read files are in standard format, they will need to be
               split into parts for bowtie to read them
  -bowmap    : Instead of a read file, a SAM mapping is given
  -bowmapDir : Instead of a read file, a directory of SAM mappings is given
  
Please see the user manual for more information about the program and options.
"""

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print __usage__
        sys.exit(0)
        
    if sys.argv[1] == '--help' or sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--h':
        print __fullUsage__
        sys.exit(0)
        
    readDir = False
    split = False
    bowmap = False
    bowmapDir = False
    buildSam = False
        
    for i in range(1,len(sys.argv)):
        if sys.argv[i][0] == '-':
            if sys.argv[i][1:] == 'readDir':
                readDir = True
                bowmap = False
                bowmapDir = False
            elif sys.argv[i][1:] == 'bowmap':
                readDir = False
                bowmap = True
                bowmapDir = False
            elif sys.argv[i][1:] == 'bowmapDir':
                readDir = False
                bowmap = False
                bowmapDir = True
            elif sys.argv[i][1:] == 'split':
                split = True
            else:
                print 'Could not recognize option : ' + sys.argv[i]
                print '\nUsage:\n'
                print __usage__
                sys.exit(1)
                
    assemFile = sys.argv[-3]
    input2 = sys.argv[-2]
    outdir = sys.argv[-1]
    
    subprocess.call("make")
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(outdir + "/indexes"):
        os.mkdir(outdir + "/indexes")

    if readDir == True:
        print "readDir not implemented yet"
        sys.exit(0)
    elif bowmap == True:
        print "bowmap not implemented yet"
        sys.exit(0)
    elif bowmapDir == True:
        print "bowmapDir not implemented yet"
        sys.exit(0)
    else:
        
        t0 = time.time()
        
        if split == True:
            print commands.getoutput("./readFileSplitter " + input2)
            commands.getoutput("mv part1_" + input2 + " " + outdir)
            commands.getoutput("mv part2_" + input2 + " " + outdir)
            print "Split took", time.time() - t0
            t0 = time.time()
        else:
            print "no split not implemented yet"
            sys.exit(0)
            
        print commands.getoutput("bowtie-build " + assemFile + " " + assemFile.split('/')[-1])
        print commands.getoutput("mv *.ebwt " + outdir + "/indexes")
        
        print "bowtie-build took", time.time() - t0
        t0 = time.time()
        
        os.chdir(outdir)
            
        print commands.getoutput("bowtie -t -I 0 -X 1000 --rf -a -l 10 -v 3 -S --sam-nohead " + assemFile.split('/')[-1] + " -1 part1_" + input2 + " -2 part2_" + input2 + " e_coli_bow_rf.map")
        print commands.getoutput("bowtie -t -I 0 -X 5000 --fr -a -l 10 -v 3 -S --sam-nohead " + assemFile.split('/')[-1] + " -1 part1_" + input2 + " -2 part2_" + input2 + " e_coli_bow_fr.map")
        
        print "bowtie took", time.time() - t0
        t0 = time.time()
        
        print commands.getoutput("cat e_coli_bow_rf.map e_coli_bow_fr.map > superMap.map")
        
        print "cat took", time.time() - t0
        t0 = time.time()
        
        os.chdir("..")
        
        print "./ALE " + outdir + "/superMap.map " + assemFile + " outputmake"
        print commands.getoutput("./ALE " + outdir + "/superMap.map " + assemFile + " outputmake")
        
        print "ALE took", time.time() - t0
        t0 = time.time()
        
        print commands.getoutput("mv outputmake " + outdir)
        
        print "Everything is done running! The output is in " + outdir + "/outputmake"
        print "Boot up a python instance and run plotter.py to start playing with the data!"