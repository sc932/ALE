import sys, os
import subprocess
import commands
import time

# parameter file for ALE server
# A reference genome is needed, there are 2 types
# unbuilt - a raw fasta reference, a needs to be run through bowtie-build (eliminates the possibility of mapped reads, all reads must be mapped to this new reference)
# built - a fasta file the corresponds to the file used to build the bowtie index used for mapping some (or all) of the reads
# There are 2 types of inputs for reads
# raw reads - these need to be split into parts, fed through bowtie and ran through ALE against a reference
# mapped reads - these need to be run through ALE against the reference
# Once the reads are mapped to the reference we call ale on each individually
# We will need to know thier insert distribution (mean, variance)
# Once all of the reads have run through ALE we need to recompile all of the outputs into a single output
# This output can now be read by the python plotter program.
#
# TAGS
# 
# [P]aired/[S]ingle (P implies the two fields after the read file(s) are mean and then variance of the insert)
# [I]nward/[O]utward (O implies the two fields after the read file(s) are mean and then variance of the outward insert, then the next two fields are mean and then variance of the inward insert)
# [1] file/[2] files
# [B]uilt/[U]nbuilt
# [R]aw/[M]apped
# I[L]lumina Scoring/S[A]M Scoring

#a typical file will look like
#----------------------------
#U AssemblyFile.fna
#PI1RL ReadFile.fna 300.0 10.0
#PI2RL part1_readFile.fna part2_readFile.fna 300.0 10.0
#PO1MS mappedReads.sam 3000.0 100.0 300.0 10.0


__version__ = '0.01'
__versionDate__ = '15 Feb 2011'
__author__ = 'Scott Clark <sc932 at cornell dot edu>'
__description__ = 'an alignment likelihood engine'
__usage__ = """python ALEserver.py [options] <ALE_server_param_file> <output dir>

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

class inputSet():
    def __init__(self, line):
        self.setOfMaps = []
        self.line = line
        self.info = line.split(' ')
        self.name = self.info[0]
        self.errors = False
        if self.name == 'U' or self.name == 'B': # this is an assembly line
            self.isReads = False
            self.file1 = self.info[1]
            if self.name == 'U':
                self.isBuilt = False
            else:
                self.isBuilt = True
        else:
            self.isReads = True
            if 'P' in self.name: # paired reads
                self.isPaired = True
                self.file1 = self.info[1]
                if '2' in self.name: # there are two files
                    self.file2 = self.info[2]
                    self.twoFiles = True
                elif '1' in self.name:
                    self.twoFiles = False
                else:
                    print "Error! Was expecting either 1 or 2 in parameter entry", self.name
                    self.errors = True
                if 'O' in self.name:
                    self.isInward = False
                    self.outwardMean = float(self.info[-4])
                    self.outwardStd = float(self.info[-3])
                    self.inwardMean = float(self.info[-2])
                    self.inwardStd = float(self.info[-1])
                elif 'I' in self.name:
                    self.isInward = True # wait for it...
                    self.inwardMean = float(self.info[-2])
                    self.inwardStd = float(self.info[-1])
                else:
                    print "Error! Was expecting either O or I in paired parameter entry", self.name
                    self.errors = True
            elif 'S' in self.name: # single reads
                self.isPaired = False
                self.file1 = self.info[1]
            else:
                print "Error! Was expecting reads to be P or S in parameter entry", self.name
                self.errors = True
            # mapped or unmapped?
            if 'M' in self.name:
                self.isMapped = True
            elif 'R':
                self.isMapped = False
            else:
                print "Error! Was expecting reads to be [R]aw or [M]apped in parameter entry", self.name
                self.errors = True
            # illumina scoring?
            if 'L' in self.name:
                self.isIllumina = True
            elif 'A' in self.name:
                self.isIllumina = False
            else:
                print "Error! Was expecting score to be I[L]lumina Scoring/S[A]M Scoring in parameter entry", self.name
                self.errors = True
    def printSpecs(self):
        print "line: ", self.line
        if self.isReads == False:
            print "Assembly File"
            print "File", self.file1
            if self.isBuilt:
                print "Built in bowtie"
            else:
                print "Still needs to be bowtie built"
        else:
            print "Read File"
            if self.isPaired:
                print "Paired Reads"
                if self.twoFiles:
                    print "Split accross two files:"
                    print "File1", self.file1
                    print "File2", self.file2
                else:
                    print "In a single file:"
                    print "File1", self.file1
                if self.isInward == False:
                    print "Outward"
                    print "Outward Mean", self.outwardMean
                    print "Outward Std", self.outwardStd
                else:
                    print "Inward"
                print "Inward Mean", self.inwardMean
                print "Inward Std", self.inwardStd
            else:
                print "Single reads"
                print "File1", self.file1
            if self.isIllumina:
                print "Illumina scoring"
            else:
                print "SAM scoring"
            if self.isMapped:
                print "Mapped"
            else:
                print "Unmapped"
            
        

def parseParameters(paramFile):
    fo = open(paramFile, 'r')
    readFiles = []
    assemFile = -1
    for line in fo:
        if line[0] != '#':
            tempFile = inputSet(line)
            if tempFile.errors == False:
                if tempFile.isReads:
                    readFiles.append(tempFile)
                else:
                    assemFile = tempFile
                tempFile.printSpecs()
            else:
                print "There were errors"
    return assemFile, readFiles
        
if __name__ == '__main__':
    if len(sys.argv) < 3:
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
        
    #for i in range(1,len(sys.argv)):
        #if sys.argv[i][0] == '-':
            #if sys.argv[i][1:] == 'readDir':
                #readDir = True
                #bowmap = False
                #bowmapDir = False
            #elif sys.argv[i][1:] == 'bowmap':
                #readDir = False
                #bowmap = True
                #bowmapDir = False
            #elif sys.argv[i][1:] == 'bowmapDir':
                #readDir = False
                #bowmap = False
                #bowmapDir = True
            #elif sys.argv[i][1:] == 'split':
                #split = True
            #else:
                #print 'Could not recognize option : ' + sys.argv[i]
                #print '\nUsage:\n'
                #print __usage__
                #sys.exit(0)
                
    stdWidth = 10.0
                
    outdir = sys.argv[-1]
    paramFile = sys.argv[-2]
    
    assemFile, readFiles = parseParameters(paramFile)
    
    subprocess.call("make")
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(outdir + "/indexes"):
        os.mkdir(outdir + "/indexes")
        
    if assemFile.isBuilt == False: # we need to build the index
        t0 = time.time()
        print "bowtie-build " + assemFile.file1 + " " + assemFile.file1.split('/')[-1]
        print commands.getoutput("bowtie-build " + assemFile.file1 + " " + assemFile.file1.split('/')[-1])
        print "mv *.ebwt " + outdir + "/indexes"
        print commands.getoutput("mv *.ebwt " + outdir + "/indexes")
        print "Building indexes took", time.time() - t0
        
    setOfOutputs = []
        
    for readFile in readFiles:
        if readFile.isPaired == False:
            print "Unpaired reads are not implemented yet"
        else: #paired reads
            if readFile.isMapped == False:
                if readFile.twoFiles == False:
                    # need to split the file
                    print commands.getoutput("./readFileSplitter " + readFile.file1)
                    commands.getoutput("mv part1_" + readFile.file1 + " " + outdir)
                    commands.getoutput("mv part2_" + readFile.file1 + " " + outdir)
                    print "Split took", time.time() - t0
                    t0 = time.time()
                    readFile.twoFiles = True
                    readFile.file1 = outdir + "/part1_" + readFile.file1
                    readFile.file2 = outdir + "/part2_" + readFile.file1
                # now we can run bowtie map
                if readFile.isInward == False:
                    # outward map
                    print "bowtie -t -I " + str(readFile.outwardMean - stdWitdh*readFile.outwardMean) + " -X " + str(readFile.outwardMean + stdWitdh*readFile.outwardMean) + " --fr -a -l 10 -v 3 -S --sam-nohead " + assemFile.file1.split('/')[-1] + " -1 part1_" + readFile.file1 + " -2 part2_" + readFile.file2 + " " + readFile.file.split('/')[-1] + "_fr.map"
                    print commands.getoutput("bowtie -t -I " + str(readFile.outwardMean - stdWitdh*readFile.outwardMean) + " -X " + str(readFile.outwardMean + stdWitdh*readFile.outwardMean) + " --fr -a -l 10 -v 3 -S --sam-nohead " + assemFile.file1.split('/')[-1] + " -1 part1_" + readFile.file1 + " -2 part2_" + readFile.file2 + " " + readFile.file.split('/')[-1] + "_fr.map")
                    readFile.setOfMaps.append(readFile.file.split('/')[-1] + "_fr.map")
                print "bowtie -t -I " + str(readFile.inwardMean - stdWitdh*readFile.inwardMean) + " -X " + str(readFile.outwardMean - stdWitdh*readFile.outwardMean) + " --rf -a -l 10 -v 3 -S --sam-nohead " + assemFile.split('/')[-1] + " -1 part1_" + readFile.file1 + " -2 part2_" + readFile.file2 + " " + readFile.file.split('/')[-1] + "_rf.map"
                print commands.getoutput("bowtie -t -I " + str(readFile.inwardMean - stdWitdh*readFile.inwardMean) + " -X " + str(readFile.outwardMean - stdWitdh*readFile.outwardMean) + " --rf -a -l 10 -v 3 -S --sam-nohead " + assemFile.split('/')[-1] + " -1 part1_" + readFile.file1 + " -2 part2_" + readFile.file2 + " " + readFile.file.split('/')[-1] + "_rf.map")
                readFile.setOfMaps.append(readFile.file.split('/')[-1] + "_rf.map")
            else:
                readFile.setOfMaps.append(readFile.file1)


    #if readDir == True:
        #print "readDir not implemented yet"
        #sys.exit(0)
    #elif bowmap == True:
        #print "bowmap not implemented yet"
        #sys.exit(0)
    #elif bowmapDir == True:
        #print "bowmapDir not implemented yet"
        #sys.exit(0)
    #else:
        
        t0 = time.time()
        
        #if split == True:
            #print commands.getoutput("./readFileSplitter " + input2)
            #commands.getoutput("mv part1_" + input2 + " " + outdir)
            #commands.getoutput("mv part2_" + input2 + " " + outdir)
            #print "Split took", time.time() - t0
            #t0 = time.time()
        #else:
            #print "no split not implemented yet"
            #sys.exit(0)
            
        #print commands.getoutput("bowtie-build " + assemFile + " " + assemFile.split('/')[-1])
        #print commands.getoutput("mv *.ebwt " + outdir + "/indexes")
        
        #print "bowtie-build took", time.time() - t0
        t0 = time.time()
        
        os.chdir(outdir)
            
        #print commands.getoutput("bowtie -t -I 0 -X 1000 --rf -a -l 10 -v 3 -S --sam-nohead " + assemFile.split('/')[-1] + " -1 part1_" + input2 + " -2 part2_" + input2 + " e_coli_bow_rf.map")
        #print commands.getoutput("bowtie -t -I 0 -X 5000 --fr -a -l 10 -v 3 -S --sam-nohead " + assemFile.split('/')[-1] + " -1 part1_" + input2 + " -2 part2_" + input2 + " e_coli_bow_fr.map")
        
        #print "bowtie took", time.time() - t0
        #t0 = time.time()
        
        #print commands.getoutput("cat e_coli_bow_rf.map e_coli_bow_fr.map > superMap.map")
        
        #print "cat took", time.time() - t0
        #t0 = time.time()
        
        os.chdir("..")
        
        print "./ALE " + outdir + "/superMap.map " + assemFile + " outputmake"
        print commands.getoutput("./ALE " + outdir + "/superMap.map " + assemFile + " outputmake")
        
        print "ALE took", time.time() - t0
        t0 = time.time()
        
        print commands.getoutput("mv outputmake " + outdir)
        
        print "Everything is done running! The output is in " + outdir + "/outputmake"
        print "Boot up a python instance and run plotter.py to start playing with the data!"