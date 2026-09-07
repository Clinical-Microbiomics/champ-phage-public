#!/usr/bin/python3

'''
Input:

Directory containing KMA output files.
A list of samples to be processed.
Output directory where the vOTU profiles and matrices will be written.

Parameters for filtering vOTUs: minimum number of reads, minimum query identity, and minimum genome coverage.
Process:

For each sample in the list, the script parses the corresponding KMA output file, filtering the vOTUs based on the provided parameters.
The filtered vOTUs and their read counts are stored in a dictionary, which is then written to a vOTU profile file.
The script also keeps track of the number of aligned and ignored reads for each sample.
After processing all samples, the script gathers all vOTU profiles and writes them into a matrix file.

Output:
vOTU profile files for each sample, containing the vOTUs and their read counts.
A matrix file containing the vOTU profiles of all samples with read counts of detected vOTUs.
'''

# Light and portable version of kmaprofiling.rs but in Python
from pathlib import Path
import csv
import argparse
import sys
class Samplestats:
    def __init__(self, alignedreads: int, ignoredreads: int, sample_id: str):
        self.alignedreads = alignedreads
        self.ignoredreads = ignoredreads
        self.sample_id = sample_id

# Parser KMA output files (.tsv) format
def parse_kma_mapping(filepath: str, sampleid: str, min_reads: int, min_cutoff: float, min_query_identity: float, read_statistics: Samplestats):
    """
    Parses the KMA mapping file and returns a dictionary of counts for each VOTU.

    Args:
        filepath (str): The path to the KMA mapping file.
        sampleid (str): The ID of the sample.
        min_reads (int): The minimum number of reads required for a vOTU to be considered.
        min_cutoff (float): The minimum coverage cutoff for a vOTU to be considered.
        min_query_identity (float): The minimum query identity for a vOTU to be considered.
        read_statistics (Samplestats): An instance of the Samplestats class to store the read statistics.

    Returns:
        dict: A dictionary containing the counts for each vOTU.

    Raises:
        Exception: If there is a formatting error in the KMA file.

    """
    alignedreads: int = 0
    ignoredreads: int = 0
    rescue_query_identity: float = 95.0
    rescue_template_cov: float = 30.0
    counts = dict()
    ignore_pattern: str = "HMR0"   # Ignore Genes from HMR-catalogue
    
    with open(filepath, 'r') as file:
        header = file.readline()    # Skip header
        linecount = 0
        for line in file:
            linecount += 1
            columns = line.split('\t')
            if len(columns) != 16:
                raise Exception(f"KMA Formatting error Sampleid: {sampleid} at line: {linecount} !!! - Consider reprofiling")  
            read_counts: int = int(columns[9])
            if read_counts < min_reads:
                ignoredreads += read_counts
                continue
            
            read_counts: int = read_counts // 2  # Adjust for PE-read - KMA counts for FE and RE together
            template_name: str = columns[0]
            
            if ignore_pattern in template_name:
                ignoredreads += read_counts
                continue
            
            votu_name = template_name.split()[0]
            query_identity: float = float(columns[5])
            template_cov: float = float(columns[3])
            
            if query_identity >= min_query_identity and template_cov >= min_cutoff:
                alignedreads += read_counts
                counts[votu_name] = read_counts
                continue
            
            if query_identity >= rescue_query_identity and template_cov >= rescue_template_cov:
                alignedreads += read_counts
                counts[votu_name] = read_counts
            else:
                ignoredreads += read_counts
    
    read_statistics.alignedreads = alignedreads
    read_statistics.ignoredreads = ignoredreads
    return counts


class KmaProfiling:
    """
    A class used to parse KMA output files and write vOTU profiles based on a sample list.

    Attributes
    ----------
    directory : Path
        The directory where the KMA output files are located.
    samples : list
        The list of samples to be processed.
    output : Path
        The directory where the output files will be written.
    min_reads : int
        The minimum number of reads required for a vOTU to be included in the profile.
    min_query_identity : float
        The minimum query identity required for a vOTU to be included in the profile.
    min_genome_coverage : float
        The minimum genome coverage required for a vOTU to be included in the profile.

    Methods
    -------
    parse_and_write_votu_kma_profiles():
        Parses KMA output files and writes vOTU profiles based on the sample list.
    write_votu_profile(filepath: Path, kma_profile: dict):
        Writes a vOTU profile to a file.
    write_matrix2file(filepath: Path, samples: list, counts: dict):
        Writes a matrix based on a dictionary containing rows of counts to a file.
    write_votu_matrices(matrixdir: Path):
        Gathers profiles into matrices and writes them to a file.
    """
    def __init__(self, directory: Path, samples: list, output: Path, min_reads: int, min_query_identity: float, min_genome_coverage: float):
        self.directory = directory
        self.samples = samples
        self.output = output
        self.min_reads = min_reads
        self.min_query_identity = min_query_identity
        self.min_genome_coverage = min_genome_coverage

    def parse_and_write_votu_kma_profiles(self):
        '''
        Parses KMA output files and writes vOTU profiles based on samplelist
        '''
        samplestatistics = []
        for sampleid in self.samples:
            kma_profile_path = self.directory.joinpath(f"kma.{sampleid}.tsv")
            
            if not kma_profile_path.exists:
                raise FileNotFoundError(f"KMA file {kma_profile_path} not found")
            
            sys.stdout.write(f"Parsing the file {kma_profile_path}\n")
            sstats = Samplestats(alignedreads=0, ignoredreads=0, sample_id=sampleid)

            try:
                kmaprofile = parse_kma_mapping(kma_profile_path, sampleid, self.min_reads, self.min_genome_coverage, self.min_query_identity, sstats)
            except Exception:
                raise FileNotFoundError("Error reading a KMA file")

            votu_profile_path = self.output.joinpath(f"{sampleid}.votu.profile.tsv")
            KmaProfiling.write_votu_profile(votu_profile_path, kmaprofile)
            samplestatistics.append(sstats)

        sys.stdout.write("Finished parsing and writing vOTU profiles\n")
        return samplestatistics
    
    @staticmethod
    def write_votu_profile(filepath: Path, kma_profile: dict):
        """
        read_counts = total # mapped reads to vOTU
        """
        with open(filepath, 'w') as writer:
            header = ["votu", "read_counts"]
            writer.write("\t".join(header) + "\n")
            for votu_name, counts in kma_profile.items():
                writer.write(f"{votu_name}\t{counts}\n")
        return None
    
    @staticmethod
    def write_matrix2file(filepath: Path, samples: list, counts: dict):
        """Open and write matrix based on dictionary containing rows of counts"""
        with open(filepath, 'w') as outfile:
            header = ['votu'] + samples
            outfile.write('\t'.join(header)+'\n')  # Write header
            for votu_name, votu_counts in counts.items():
                string_counts = '\t'.join(str(f) for f in votu_counts)
                outfile.write(votu_name+'\t'+string_counts+'\n')
        return None

    def write_votu_matrices(self, matrixdir: Path):
        """Gather Profiles into Matrices write to file"""
        samples = self.samples
        nsamples = len(samples)

        # Read count profiles
        sys.stdout.write("Gathering vOTU read-count profiles\n")
        votu_arrays = {}
        for sidx, sampleid in enumerate(samples):
            votu_profile_path = self.output.joinpath(f"{sampleid}.votu.profile.tsv")
            with open(votu_profile_path, "r") as votufile:
                reader = csv.reader(votufile, delimiter="\t")
                next(reader)  # Skip Header
                for line in reader:
                    votu_name = line[0]
                    read_counts = int(line[1])
                    if votu_name in votu_arrays:
                        # Key exists, modify the value
                        votu_arrays[votu_name][sidx] = read_counts
                    else:
                        # Key doesn't exist, insert a new entry
                        counts = [0] * nsamples
                        counts[sidx] = read_counts
                        votu_arrays[votu_name] = counts

        matrix_path = matrixdir.joinpath("votu.counts.mat.tsv")
        KmaProfiling.write_matrix2file(matrix_path, samples, votu_arrays)

        # Record votus profiled
        votusprofiled = list(votu_arrays.keys())
        return votusprofiled



if __name__=='__main__':

    # Argparse with output argument
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--kmadir',type=str,help='KMA files directory')
    parser.add_argument('-o','--output',type=str,help='Output vOTU counts matrix')
    args = parser.parse_args()
    
    kma_directory = Path(args.kmadir)
    matrices_out_directory = Path(args.output)
    profiles_directory = matrices_out_directory.joinpath("votu_profiles")
    profiles_directory.mkdir(parents=True, exist_ok=True)
    matrices_out_directory.mkdir(parents=True, exist_ok=True)

    # Parse samples
    samplefile_path = Path("samplelist")
    if not samplefile_path.exists():
        raise FileNotFoundError(f"Samplelist file {samplefile_path} not found")
    
    with open(samplefile_path, "r") as samplefile:
        reader = samplefile.readlines()
        samples = [line.split("\t")[0].strip() for line in reader]

    # Initialise parameters for kma profiling
    kma_profiling = KmaProfiling(
        directory=kma_directory,
        samples=samples.copy(),
        output=profiles_directory,
        min_reads=10,
        min_query_identity=99.0,
        min_genome_coverage=5.0,
    )

    # Parse and Write votu profiles after filtering
    try:
        kmastats = kma_profiling.parse_and_write_votu_kma_profiles()
    except Exception as err:
        raise FileExistsError(f"Error parsing votu profiles {err}")

    # Write count matrix by parsing all votu profiles writing to out directory
    try:
        votusprofiled = kma_profiling.write_votu_matrices(matrices_out_directory)
    except Exception as err:
        raise Exception(f"Error when writing vOTU matrices -> {err}")
