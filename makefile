.DEFAULT_GOAL := model

# FASTA and VCF options
# --------------------------------------------------------------------------------------
SAVE = /opt/UPV/TFG/src/example/
VCF = -vcf $(SAVE)datosR1.vcf
VCF_2 = -vcf /home/alejandro/Descargas/datosvcf/RP924_9589186940.vcf
FASTA_2 = -fasta /home/alejandro/GRCh37.p13.genome.fa.gz
FASTA = -fasta $(SAVE)hg19.fa.gz
FILE_OPTIONS = $(VCF) $(FASTA) -s $(SAVE)

# Sequences options
# --------------------------------------------------------------------------------------
LENGTH_PREFIX = 20
LENGTH_SUFFIX = 20
SEQUENCE_OPTIONS = -p_p $(LENGTH_PREFIX) -p_s $(LENGTH_SUFFIX)

# PARSER options
# --------------------------------------------------------------------------------------
PARSER_TYPE = m
PARSER_ADD_ORIGINAL = -ao
PARSER_OPTIONS = -p $(PARSER_TYPE) $(PARSER_ADD_ORIGINAL)

# Model options
# --------------------------------------------------------------------------------------
K = 3
MODEL = ktss
KTSS_OPTIONS = -m $(MODEL) -k $(K)

# VALIDATOR options
# --------------------------------------------------------------------------------------
ADD_MUTATION_VALIDATOR = -amval
ADD_ORIGNAL_VALIDATOR = -aoval
SAVE_DISTANCES = -sd
VALIDATOR_OPTIONS = $(ADD_ORIGNAL_VALIDATOR) -min $(ADD_MUTATION_VALIDATOR) $(SAVE_DISTANCES)

# General model options
# --------------------------------------------------------------------------------------
RATIO = 0.90
STEPS = 1
INIT = python3 init.py
TEST =
GENERAL_OPTIONS = $(INIT) -r $(RATIO) -steps $(STEPS) -o pm $(FILE_OPTIONS) $(TEST)

build-docs:
	pdoc --html src --force
build-docs-pdf:
	pdoc --pdf src --force
test:
	poetry run pytest
black:
	poetry run black .
parser:
	$(INIT) $(PARSER_OPTIONS) $(FILE_OPTIONS)
model:
	$(GENERAL_OPTIONS) $(SEQUENCE_OPTIONS) $(PARSER_OPTIONS) $(KTSS_OPTIONS) $(VALIDATOR_OPTIONS)
