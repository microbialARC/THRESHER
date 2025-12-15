import re
# Replace any non-alphanumeric characters (except dash/underscore) with underscore.
def parse_genome_name(original_name: str) -> str:
    parsed_name = re.sub(r'[^a-zA-Z0-9._-]', '_', original_name)
    # Collapse multiple underscores
    parsed_name = re.sub(r'_+', '_', parsed_name)
    # If the sample name is "Reference", change it to "Reference_Genome" to avoid conflict in snippy.
    if parsed_name == "Reference":
        parsed_name = "Reference_Genome"
    return parsed_name.strip('_')
