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

    def genotype(self, extracted_fq, output_dir, genes, threads=8, verbose=True, single=False):
        """
        Run the 'genotype' command of arcasHLA.

        Parameters:
        - extracted_fq (str): Path to the extracted FASTQ file.
        - output_dir (str): Directory to save the genotyping results.
        - genes (str): Comma-separated list of genes to include.
        - threads (int): Number of threads to use.
        - verbose (bool): Whether to enable verbose output.
        - single (bool): Whether to enable single-end mode.

        Returns:
        - None
        """
        verbose_flag = "-v" if verbose else ""
        single_flag = "--single" if single else ""
        command = f"{self.arcas_hla_path} genotype {extracted_fq} -o {output_dir} -g {genes} -t {threads} {verbose_flag} {single_flag}"
        self.run_command(command)

    def partial_genotype(self, extracted_fq, genes, genotype_file, output_dir, threads=8, verbose=True, single=True):
        """
        Run the 'partial' command of arcasHLA.

        Parameters:
        - extracted_fq (str): Path to the extracted FASTQ file.
        - genes (str): Comma-separated list of genes to include.
        - genotype_file (str): Path to the genotype JSON file.
        - output_dir (str): Directory to save the partial genotyping results.
        - threads (int): Number of threads to use.
        - verbose (bool): Whether to enable verbose output.
        - single (bool): Whether to enable single-end mode.

        Returns:
        - None
        """
        verbose_flag = "-v" if verbose else ""
        single_flag = "--single" if single else ""
        command = (f"{self.arcas_hla_path} partial {extracted_fq} -g {genes} -G {genotype_file} "
                   f"-o {output_dir} -t {threads} {verbose_flag} {single_flag}")
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
