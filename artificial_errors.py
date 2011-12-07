#!/usr/bin/python

# (C) 2011 Scott Clark

__usage__ = """Usage: ./artificial_errors.py [-options] <inputfile.fna>

where basic options are:
  -h      : show brief help on version and full usage
"""
__full_usage__="""Usage: ./artificial_errors.py [-options] <inputfile.fna>

where basic options are:
  -h      : show brief help on version and full usage

parameter options accepting <i>ntegers and <s>trings (default):
  Note: transformations will be made left to right
  -ase <i> <i> : add substitution error at <location> for <length>
  -ade <i> <i> : add deletion error at <location> for <length>
  -aie <i> <i> : add insertion error at <location> for <length>
  -inv <i> <i> : add inversion error at <location> for <length>
  -cip <i> <i> : copy part of the assembly at <location> for <length>
  -trp <i>     : transpose assembly around <pivot>
  -ab  <i>     : add a break (split into 2 contigs) at <location>
  -o   <s>     : output file name (error_ + inputfile.fna)
"""

import sys

inversion_error = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

rotate_error = {'A':'T', 'T':'C', 'C':'G', 'G':'A'}

def add_sub_error(assembly, location, length):
    """add a substitution error to <assembly> at <location> for a set <length>"""
    for i in range(length):
        assembly[location + i] = rotate_error[assembly[location + i]]

def add_insertion_error(assembly, location, length):
    """add a insertion error to <assembly> at <location> for a set <length>"""
    for _ in range(length):
        assembly.insert(location, 'A')

def add_deletion_error(assembly, location, length):
    """add a deletion error to <assembly> at <location> for a set <length>"""
    for _ in range(length):
        assembly.pop(location)

def transpose_assembly(assembly, pivot):
    """transpose <assembly> about some <pivot>"""
    temp = assembly[:pivot][:]
    assembly[:len(assembly)-pivot] = assembly[-len(assembly)+pivot:][:]
    assembly[-pivot:] = temp[:]

def add_inversion_error(assembly, location, length):
    """add a inversion error to <assembly> at <location> for a set <length>"""
    for i in range(length):
        assembly[location + i] = inversion_error[assembly[location + i]]

def read_in_assembly(assembly_file):
    """read a fasta file <assembly_file> into a list of bases and return it"""
    input_file = open(assembly_file, 'r')
    assembly = []
    for line in input_file:
        if line[0] != '>': # ignore comments
            assembly.extend(list(line[:-1])) # add it to the assembly, drop the \n
    input_file.close()
    return assembly

def copy_in_place(assembly, location, length):
    """copy a section of <assembly> at <location> for a given <length>"""
    end = assembly[location:][:]
    start = assembly[:location][:]
    copy = assembly[location:location+length][:]
    assembly = start
    assembly.extend(copy)
    assembly.extend(end)
    return assembly

def add_break(assembly, location):
    """add a break in <assembly> at <location>"""
    assembly.insert(location, 'break')

def output_assembly(file_name, assembly):
    """output a list of bases as an assembly file <file_name>"""
    output_file = open(file_name, 'w')
    contig = 0
    output_file.write('>' + file_name + ' ' + str(contig) + '\n')
    output_line = ""
    base_on = 0
    for base in assembly:     
        if base == 'break':
            if output_line != "":
                output_file.write(output_line + '\n')
                output_line = ""
            contig += 1
            output_file.write('>' + file_name + ' ' + str(contig) + '\n')
            base_on = 0
        else:
            base_on += 1
            output_line += base
            if not base_on%70:
                output_line += '\n'
                output_file.write(output_line)
                output_line = ""
    if output_line != "":
        output_file.write(output_line + '\n')
    output_file.close()

def main():
    """read in an assembly file and transform it based on the command line options
       then output it again in another fasta file, see __full_usage__"""
    if len(sys.argv) < 2:
        print __usage__
        sys.exit(0)
        
    if sys.argv[1] == '--help' or sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--h':
        print __full_usage__
        sys.exit(0)

    assembly = read_in_assembly(sys.argv[-1])

    output_file = 'errors_' + sys.argv[-1]

    arg_on = 1
    while(arg_on + 1 < len(sys.argv)):            
        if sys.argv[arg_on] == '-ase':
            location = int(sys.argv[arg_on + 1])
            length = int(sys.argv[arg_on + 2])
            add_sub_error(assembly, location, length)
            arg_on += 3
        elif sys.argv[arg_on] == '-aie':
            location = int(sys.argv[arg_on + 1])
            length = int(sys.argv[arg_on + 2])
            add_insertion_error(assembly, location, length)
            arg_on += 3
        elif sys.argv[arg_on] == '-ade':
            location = int(sys.argv[arg_on + 1])
            length = int(sys.argv[arg_on + 2])
            add_deletion_error(assembly, location, length)
            arg_on += 3
        elif sys.argv[arg_on] == '-inv':
            location = int(sys.argv[arg_on + 1])
            length = int(sys.argv[arg_on + 2])
            add_inversion_error(assembly, location, length)
            arg_on += 3
        elif sys.argv[arg_on] == '-cip':
            location = int(sys.argv[arg_on + 1])
            length = int(sys.argv[arg_on + 2])
            assembly = copy_in_place(assembly, location, length)
            arg_on += 3
        elif sys.argv[arg_on] == '-trp':
            pivot = int(sys.argv[arg_on + 1])
            transpose_assembly(assembly, pivot)
            arg_on += 2
        elif sys.argv[arg_on] == '-ab':
            location = int(sys.argv[arg_on + 1])
            add_break(assembly, location)
            arg_on += 2
        elif sys.argv[arg_on] == '-o':
            output_file = sys.argv[arg_on + 1]
            arg_on += 2
        else:
            print "Did not recognize command line argument %s." % sys.argv[arg_on]
            print "Try -h for help."
            exit(0)

    output_assembly(output_file, assembly)
    print "Output assembly in %s" % output_file

if __name__ == '__main__':
    main()

