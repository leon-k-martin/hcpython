import re
from pathlib import Path

# Load the Snakefile content
snakefile_path = Path("/Users/leonmartin_bih/projects/hcpython/Snakefile")
snakefile_content = snakefile_path.read_text()

# Regular expressions to find inputs and outputs in the Snakefile
input_pattern = re.compile(r"input:\s*(.+)", re.MULTILINE)
output_pattern = re.compile(r"output:\s*(.+)", re.MULTILINE)
file_pattern = re.compile(r'\"([^\"]+)\"')

# Function to extract file paths from a given match object
def extract_files(block):
    files = []
    matches = file_pattern.findall(block)
    files.extend(matches)
    return files

# Extract all inputs and outputs
inputs = re.findall(input_pattern, snakefile_content)
outputs = re.findall(output_pattern, snakefile_content)

# Flatten the list of inputs and outputs
input_files = [item for sublist in inputs for item in extract_files(sublist)]
output_files = [item for sublist in outputs for item in extract_files(sublist)]

# Convert lists to sets for easier comparison
input_files_set = set(input_files)
output_files_set = set(output_files)

# Determine input files not created by any rule
uncreated_files = input_files_set - output_files_set

# Print the result
print("Input files not created by any rule:")
for file in sorted(uncreated_files):
    print(file)