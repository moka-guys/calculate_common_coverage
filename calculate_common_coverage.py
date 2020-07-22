import subprocess
import argparse
import sys

def cli_arguments(args):
    """Parses command line arguments.
    Args:
        args: A list containing the expected commandline arguments. Example:
            ['calculate_common_coverage.py', 
            '-r', '20', 
            '-b', '/path/to/bedfile.bed',
            '-m',
            '-q', '20',
            '-o', '/path/to/outputfile',
            '-t', '/path/to/sample1_markdup_sorted.bam',
            '-u', '/path/to/sample2_markdup_sorted.bam', 
            '-v', '/path/to/sample3_markdup_sorted.bam']
    Returns:
        An argparse.parser object with methods named after long-option command-line arguments. Example:
            --runfolder "media/data1/share/runfolder" --> parser.parse_args(args).runfolder
    """
    # Define arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--readdepth', required=True, default=20, help='The required read depth, default = 20')
    parser.add_argument('-b', '--bedfile', required=True, help='path to bedfile, required')
    parser.add_argument('-t', '--bamfile1', required=True, help='path to bamfile1, required')
    parser.add_argument('-u', '--bamfile2', required=True, help='path to bamfile2, required')
    parser.add_argument('-v', '--bamfile3', required=False, default=" ", help='path to bamfile3 - optional')
    parser.add_argument('-q', '--minbasequal', required=False, default=" ", help='minimum base quality to be counted - optional')
    parser.add_argument('-m', '--countoverlapreadsonce', action='store_true', default=False, help='-m does not count overlapping mate reads more than once')
    parser.add_argument('-o', '--output_file', required=True, help='path to output file, required')
    
    # Collect arguments and return
    return parser.parse_args(args)

class trio_coverage():
    """
    This class takes a BED file, two or three BAM files and arguments for sambamba coverage.
    It will calculate the read depth at each base for each sample and report if this base is covered sufficiently in all samples
    The read depth is calculated using sambamba, mirroring the existing coverage calculation.
    The % of each gene covered sufficiently is calculated by counting the number of bases covered sufficiently and those not covered sufficiently.
    A test is performed to assess if any amplicons in the gene overlap, as overlapping regions would result in bases being counted twice, skewing the % of the geen which is covered sufficiently
    The result is output in the same format as the chanjo_txt files, used to import WES coverage into ngscoverage table
    """
    def __init__(self, min_coverage, bedfile_path, bamfile1_path, bamfile2_path, bamfile3_path, output_file_path, countoverlapreadsonce,minbasequal):
        self.min_coverage = min_coverage
        self.countoverlapreadsonce = countoverlapreadsonce
        self.minbasequal = minbasequal
        self.coverage_dict = {}
        self.bamfile1_path = bamfile1_path
        self.bamfile2_path = bamfile2_path
        self.bamfile3_path = bamfile3_path
        self.bedfile_path = bedfile_path
        self.output_file_path = output_file_path

        print(
            "min coverage = {}\n".format(self.min_coverage),
            "count overlap reads once = {}\n".format(self.countoverlapreadsonce),
            "minbasequal = {}\n".format(self.minbasequal),
            "bamfile1_path = {}\n".format(self.bamfile1_path),
            "bamfile2_path = {}\n".format(self.bamfile2_path),
            "bamfile3_path = {}\n".format(self.bamfile3_path),
            "bedfile = {}\n".format(self.bedfile_path),
            "output path = {}".format(self.output_file_path)
        )
    
    def get_amplicons(self):
        """
        Parse the BED file to populate a dictionary with an entry for each entrez gene id and the value a list of amplicons
        Sort the amplicons within each gene on chromsome and start to ensure any overlaps would be identified
        """
        # create a local dict, before sorting
        coverage_dict = {}
        with open(self.bedfile_path,'r') as bedfile:
            for line in bedfile.readlines():
                # create dict with entrez_gene_id as key and list of tuples describing amplicon regions
                chr, start, stop, fourth, fifth, sixth, seventh, entrez = line.rstrip().split("\t")
                coverage_dict.setdefault(entrez, []).append((chr, start, stop))
        for entrez in coverage_dict:
            # for each entrezid sort the list on chrom and start and append sorted list to the class-wide dictionary
            self.coverage_dict[entrez] = sorted(coverage_dict[entrez], key=lambda x: (x[0],x[1]))

    def build_sambamba_opts(self):
        """
        Build and return a string of arguments for the sambamba command.
        Uses countoverlapreadsonce and minbasequal arguments provided as inputs to the script (and captured as class-wide variables)
        """
        opts_string = ""
        # if boolean arg provided add -m flag
        if self.countoverlapreadsonce:
            opts_string += " -m "
        # if minbasequal not given to script given default value of " "
        # if this is provided  add to input in format -q=10
        if self.minbasequal != " ":
            opts_string += " -q=" + self.minbasequal
        return opts_string

    def run_sambamba(self, chr, start, stop, bamfile1, bamfile2, bamfile3, opts_list):
        """
        Build and run the sambamba command for up to three BAM files for one amplicon, provided in format of chr:start-stop
        Sambamba depth base command is used to output per base coverage
        Output contains a mix of warning/error messages and base level coverage.
        Sambamba outputs to stdout which is split into a list of lines.
        Each line is split on tab, created a list for each line 
        Inputs:
            chr - chromosome of amplicon
            start - start coord of amplicon
            stop - stop coord of amplicon
            bamfile1 - path to bamfile1
            bamfile2 - path to bamfile2
            bamfile3 - path to bamfile3 (if not provided to script default value is " ")
            opts_list - string of options created by build_sambamba_opts()
        Returns:
            List of lists (list where each item in list is a line, and each element is itself a list, splitting the line on tabs 
            eg [chrom, pos, cov, a, c, g, t, deleted, refskip, sample])
        """
        ## sambamba settings
        # -c minimum mean coverage for output (default: 0 for region/window, 1 for base)
        # -L single region in form chr:beg-end
        # -m detect overlaps of mate reads and handle them on per-base basis (if given)
        # -q only count bases with a base quality score greater or equal to this value (if given)
        # -F filter reads using the filter eg 'mapping_quality >= 20'
        sambamba_command = "'/resources/sambamba-0.7.1-linux-static'" \
            + " depth base {} -c 0 -L {}:{}-{} -F 'mapping_quality >= 20' {} {} {} \
            ".format(opts_list, chr, start, stop, bamfile1, bamfile2, bamfile3)
        sambamba_out = subprocess.run(sambamba_command, shell=True, capture_output=True, universal_newlines=True)
        sambamba_list = []
        # split stdout into lines (using newlines), and each line on tab
        for i in sambamba_out.stdout.split("\n"):
            sambamba_list.append(i.split("\t"))
        return sambamba_list
        
    
    def test_amplicons_in_gene(self,entrez):
        """
        Look for overlap between adjacent regions within a gene
        The entrez gene id is the key to self.coverage_dict, which has a value of sorted list of tuples of amplicon coords
        This function parses these and assesses for overlap between adjacent regions
        Input: 
            Entrez gene id
        Returns:
            Error is raised if overlapping amplicons are found 
        """
        # loop through list of amplicons
        for amplicon in range(0, len(self.coverage_dict[entrez])):
            chr1, start1, stop1 = self.coverage_dict[entrez][amplicon]
            # if not the last amplicon
            if amplicon < len(self.coverage_dict[entrez]) - 1:
                # find the coordinates of the next amplicon in bed
                chr2, start2, stop2 = self.coverage_dict[entrez][amplicon + 1]
                # check the amplicons don't overlap
                if chr1 == chr2 and int(start2) < int(stop1):
                    raise AssertionError("BED file has overlapping regions within an entrez gene id")

    def parse_sambamba_output(self, sambamba_list):
        """
        Parses the sambamba output for an amplicon (returned by run_sambamba()), extracting only the relevant lines and fields.
        Line of interest contains the read depth at one base for one sample. Other lines include warning messages.
        If it's a line containing a read depth of interest, extract the read depth, and coordinate 
        Inputs:
            A list of lists produced by run_sambamba()
            Each item in the list is a line output from sambamba.
            Each line has been split on tab into another list.
            Lines of interest are lists with following format [chrom, pos, cov, a, c, g, t, deleted, refskip, sample]
        Returns:
            An amplicon dictionary with the pos as a key (in form chr:pos) and the value a list of read depths at that position
        """
        amplicon_dict = {}
        for line in sambamba_list:
            # 'skip header / empty lines
            if not line[0] == "REF" and len(line) == 10:
                # assign list to variables
                chrom, pos, cov, a, c, g, t, deleted, refskip, sample = line
                # build coords
                coord = chrom + ":" + pos
                # add to dictionary, with the coord as the key and a list of coverages as value
                amplicon_dict.setdefault(coord,[]).append(int(cov))
        return amplicon_dict

    
    def calculate_gene_coverage(self):
        """
        This function loops through self.coverage_dict, which describes the regions in the BED file.
        For each gene the number of bases which are, and are not covered sufficiently in all samples are counted.
        This is done by calling sambamba on each exon/amplicon at a time
        Checks are performed to make sure we are taking into account all bases.
        The entrezgeneID and % covered in all samples is written to the output file provided as an argument to the script.
        """
        with open(self.output_file_path,'w') as output_file:
            # for each gene
            for entrez in self.coverage_dict:
                self.test_amplicons_in_gene(entrez)
                # set up counts for bases which are or are not covered sufficiently in ALL BAMS
                gene_level_covered_pass_count = 0
                gene_level_covered_fail_count = 0
                # for each amplicon in BED file
                for amplicon in self.coverage_dict[entrez]:
                    # create counts for each amplicon
                    amplicon_base_ok = 0
                    amplicon_base_fail = 0    
                    # unpack coordinate tuple
                    chr,start,stop = amplicon
                    # set length of amplicon, taking into account the zero based, half open BED file 
                    amplicon_length = int(stop)-int(start)+1
                    # build_sambamba_opts() returns list of sambamba arguments
                    # pass sambamba arguments, coordinates and BAM files to run_sambamba() 
                    # pass output of run_sambamba() into parse_sambamba_output() which returns a dictionary of read depths at each base in the amplicon
                    amplicon_dict = self.parse_sambamba_output(
                        self.run_sambamba(chr, start, stop, self.bamfile1_path, self.bamfile2_path, self.bamfile3_path, self.build_sambamba_opts())
                        )
                    
                    # pass this dictionary determining if all samples are covered sufficiently and add to count 
                    for base in amplicon_dict:
                        # if lowest covered base at this base is less than minimum coverage level add to relevant count.
                        if min(amplicon_dict[base]) < int(self.min_coverage):
                            amplicon_base_fail +=1
                        else:
                            amplicon_base_ok +=1
                        
                    # check expected number of bases in the amplicon have been checked
                    assert amplicon_base_ok + amplicon_base_fail == amplicon_length
                    # add amplicon counts to gene wide counts
                    gene_level_covered_pass_count += amplicon_base_ok
                    gene_level_covered_fail_count += amplicon_base_fail
                # calculate percentage of gene covered sufficiently (pass count/total base count)*100
                gene_percent = str(gene_level_covered_pass_count / (gene_level_covered_pass_count + gene_level_covered_fail_count) * 100)
                # write gene level results to output file
                output_file.write("\t".join([entrez, gene_percent, "NULL\n"]))
        
    

def main(args):
    # Get command line arguments
    parsed_args = cli_arguments(args)
    tc = trio_coverage(min_coverage=parsed_args.readdepth, \
        bedfile_path=parsed_args.bedfile, \
        bamfile1_path=parsed_args.bamfile1, \
        bamfile2_path=parsed_args.bamfile2, \
        bamfile3_path=parsed_args.bamfile3, \
        output_file_path=parsed_args.output_file, \
        countoverlapreadsonce=parsed_args.countoverlapreadsonce, \
        minbasequal=parsed_args.minbasequal)
    tc.get_amplicons()
    tc.calculate_gene_coverage()

if __name__ == '__main__':
    main(sys.argv[1:])

    #python3 /home/aled/Documents/DNAnexus_apps/dnanexus_calculate_common_coverage/resources/home/dnanexus/calculate_common_coverage.py -r 20 -b /home/aled/Documents/mokabed/LiveBedfiles/Pan1449dataSambamba.bed -t /home/aled/Documents/200416_trio_coverage/NGS330_34_NA24385_CN_M_Exome_Pan3174_S2_markdup_sorted.bam -u /home/aled/Documents/200416_trio_coverage/NGS330_35_NA24149_CN_M_Exome_Pan3174_S3_markdup_sorted.bam -v /home/aled/Documents/200416_trio_coverage/NGS330_37_222577_JC_M_WES55_Pan2835_S5_markdup_sorted.bam -m -q 0 -o ~/test_output
    