.DEFAULT_GOAL := model

# FASTA and VCF options
# --------------------------------------------------------------------------------------
SAVE = /opt/UPV/TFG/src/example/
VCF = -vcf $(SAVE)datosR1.vcf
FASTA = -fasta $(SAVE)hg19.fa.gz
FILE_OPTIONS = $(VCF) $(FASTA) -s $(SAVE)

# Sequences options
# --------------------------------------------------------------------------------------
LENGTH = 20
SEQUENCE_OPTIONS = -p_p $(LENGTH) -p_s $(LENGTH)

# PARSER options
# --------------------------------------------------------------------------------------
PARSER_TYPE = m
PARSER_ADD_ORIGINAL = -ao
PARSER_OPTIONS = -p $(PARSER_TYPE) $(PARSER_ADD_ORIGINAL) -amto

# Model options
# --------------------------------------------------------------------------------------
K = 3
MODEL = ktss
KTSS_OPTIONS = -m $(MODEL) -k $(K)

# VALIDATOR options
# --------------------------------------------------------------------------------------
ADD_MUTATION_VALIDATOR = -amval
ADD_ORIGNAL_VALIDATOR = -aoval
VALIDATOR_OPTIONS = $(ADD_ORIGNAL_VALIDATOR) -min $(ADD_MUTATION_VALIDATOR)

# General model options
# --------------------------------------------------------------------------------------
RATIO = 0.90
STEPS = 1
INIT = python3 init.py
GENERAL_OPTIONS = $(INIT) -r $(RATIO) -steps $(STEPS) -o pm $(FILE_OPTIONS)

build-docs:
	pdoc --html --config show_source_code=False src --force
build-docs-pdf:
	pdoc --pdf --config show_source_code=False src --force
test:
	poetry run pytest
black:
	poetry run black .
parser:
	$(INIT) $(PARSER_OPTIONS) $(FILE_OPTIONS)
model:
	$(GENERAL_OPTIONS) $(SEQUENCE_OPTIONS) $(PARSER_OPTIONS) $(KTSS_OPTIONS) $(VALIDATOR_OPTIONS)
