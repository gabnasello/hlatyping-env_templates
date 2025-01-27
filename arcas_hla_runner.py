import subprocess
import os

class ArcasHLARunner:
    def __init__(self, conda_env="arcas-hla", arcas_hla_path="arcasHLA/arcasHLA"):
        """
        Initialize the ArcasHLARunner class.

        Parameters:
        - conda_env (str): Name of the Conda environment to use.
        - arcas_hla_path (str): Path to the arcasHLA executable.
        """
        self.conda_env = conda_env
        self.arcas_hla_path = arcas_hla_path

    def set_conda_env(self, env_name):
        """Update the Conda environment to use."""
        self.conda_env = env_name

    def set_arcas_hla_path(self, path):
        """Update the path to the arcasHLA executable."""
        self.arcas_hla_path = path
        version_file = os.path.join(os.path.dirname(path), "dat/IMGTHLA/release_version.txt")
        if os.path.exists(version_file):
            with open(version_file, "r") as f:
                print(f"Current reference database version: {f.read().strip()}")
        else:
            print("Error: Reference database version file not found. Please update the database before proceeding.")

    def run_command(self, command):
        """
        Run a terminal command within the specified Conda environment and print its output live.
    
        Parameters:
        - command (str): The command to execute.
    
        Returns:
        - None
        """
        full_command = f"conda run -n {self.conda_env} {command}"
        print(f"Executing command: {full_command}")  # Print the command before execution
        process = subprocess.Popen(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(process.stdout.readline, b""):
            print(line.decode('utf-8'), end="")
        process.stdout.close()
        process.wait()


    def extract_reads(self, bam_file, output_dir, threads=8, verbose=True, single=False):
        """
        Run the 'extract' command of arcasHLA.

        Parameters:
        - bam_file (str): Path to the input BAM file.
        - output_dir (str): Directory to save the extracted reads.
        - threads (int): Number of threads to use.
        - verbose (bool): Whether to enable verbose output.
        - single (bool): Whether to enable single-end mode.

        Returns:
        - None
        """
        verbose_flag = "-v" if verbose else ""
        single_flag = "--single" if single else ""
        command = f"{self.arcas_hla_path} extract {bam_file} -o {output_dir} -t {threads} {verbose_flag} {single_flag}"
        self.run_command(command)

    def genotype(
        self,
        extracted_fq,
        output_dir,
        genes,
        population="prior",
        min_count=75,
        tolerance=1e-7,
        max_iterations=1000,
        drop_iterations=None,
        drop_threshold=0.1,
        zygosity_threshold=0.15,
        log_file=None,
        temp_dir="/tmp",
        keep_files=False,
        threads=1,
        verbose=False,
        single=False,
        avg_fragment_length=200,
        std_fragment_length=20
    ):
        """
        Run the 'genotype' command of arcasHLA.
    
        Parameters:
        - extracted_fq (str): Path to the extracted FASTQ file.
        - output_dir (str): Directory to save the genotyping results.
        - genes (str): Comma-separated list of genes to include.
        - population (str): Sample population (e.g., "asian_pacific_islander", "black", "caucasian", etc.).
        - min_count (int): Minimum gene read count required for genotyping.
        - tolerance (float): Convergence tolerance for transcript quantification.
        - max_iterations (int): Maximum number of iterations for transcript quantification.
        - drop_iterations (int): Number of iterations before dropping low support alleles.
        - drop_threshold (float): Proportion of maximum abundance an allele needs to not be dropped.
        - zygosity_threshold (float): Threshold for determining zygosity.
        - log_file (str): Path to log file for run summary.
        - temp_dir (str): Temporary directory to use.
        - keep_files (bool): Whether to keep intermediate files.
        - threads (int): Number of threads to use.
        - verbose (bool): Whether to enable verbose output.
        - single (bool): Whether to enable single-end mode.
        - avg_fragment_length (int): Estimated average fragment length for single-end reads.
        - std_fragment_length (int): Estimated standard deviation of fragment length.
    
        Returns:
        - None
        """
        verbose_flag = "-v" if verbose else ""
        single_flag = "--single" if single else ""
        keep_files_flag = "--keep_files" if keep_files else ""
        drop_iterations_flag = f"--drop_iterations {drop_iterations}" if drop_iterations is not None else ""
        log_flag = f"--log {log_file}" if log_file else ""
    
        command = (
            f"{self.arcas_hla_path} genotype {extracted_fq} "
            f"-o {output_dir} -g {genes} -p {population} --min_count {min_count} "
            f"--tolerance {tolerance} --max_iterations {max_iterations} {drop_iterations_flag} "
            f"--drop_threshold {drop_threshold} --zygosity_threshold {zygosity_threshold} "
            f"{log_flag} --temp {temp_dir} {keep_files_flag} -t {threads} {verbose_flag} {single_flag} "
            f"-l {avg_fragment_length} -s {std_fragment_length}"
        )
    
        self.run_command(command)


    def partial_genotype(self, extracted_fq, genotype_file, output_dir, genes, threads=8, verbose=True, single=True):
        """
        Run the 'partial' command of arcasHLA.
    
        Parameters:
        - extracted_fq (str): Path to the extracted FASTQ file.
        - genotype_file (str): Path to the genotype JSON file.
        - output_dir (str): Directory to save the partial genotyping results.
        - genes (str): Comma-separated list of genes to include.
        - threads (int): Number of threads to use.
        - verbose (bool): Whether to enable verbose output.
        - single (bool): Whether to enable single-end mode.
    
        Returns:
        - None
        """
        verbose_flag = "-v" if verbose else ""
        single_flag = "--single" if single else ""
        command = (f"{self.arcas_hla_path} partial {extracted_fq} -G {genotype_file} -o {output_dir} "
                   f"-g {genes} -t {threads} {verbose_flag} {single_flag}")
        
        self.run_command(command)        


    def update_reference(self):
        """
        Update the IPD-IMGT/HLA database to the latest version.

        Returns:
        - None
        """
        command = f"{self.arcas_hla_path} reference --update"
        self.run_command(command)
        version_file = os.path.join(os.path.dirname(self.arcas_hla_path), "dat/IMGTHLA/release_version.txt")
        if os.path.exists(version_file):
            with open(version_file, "r") as f:
                print(f"Updated reference database version: {f.read().strip()}")

    def rollback_reference(self, version):
        """
        Roll back the IPD-IMGT/HLA database to a specific version.

        Parameters:
        - version (str): The version to roll back to (e.g., "3.24.0").

        Returns:
        - None
        """
        command = f"{self.arcas_hla_path} reference --version {version}"
        self.run_command(command)
        version_file = os.path.join(os.path.dirname(self.arcas_hla_path), "dat/IMGTHLA/release_version.txt")
        if os.path.exists(version_file):
            with open(version_file, "r") as f:
                print(f"Rolled back to reference database version: {f.read().strip()}")

# Example usage
if __name__ == "__main__":
    runner = ArcasHLARunner()
    runner.set_arcas_hla_path("arcasHLA/arcasHLA")

    # Example workflow
    # Step 1: Extract reads
    # runner.extract_reads("samples/GC130443/GC130443.bam", "samples/GC130443/", threads=8, verbose=True, single=True)

    # Step 2: Genotype
    # runner.genotype("samples/GC130443/GC130443.extracted.fq.gz", "samples/GC130443/results/", genes="A,B,C,DPB1,DQB1,DQA1,DRB1", threads=8, verbose=True, single=True)

    # Step 3: Partial Genotyping
    # runner.partial_genotype("samples/GC130443/GC130443.extracted.fq.gz", "A,B,C,DPB1,DQB1,DQA1,DRB1", \
    #                         "samples/GC130443/output_genotype/GC130443.genotype.json", \
    #                         "samples/GC130443/output_partial", threads=8, verbose=True, single=True)

    # Step 4: Update reference database
    # runner.update_reference()

    # Step 5: Rollback reference database to a specific version
    # runner.rollback_reference("3.24.0")
