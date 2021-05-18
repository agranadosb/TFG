.DEFAULT_GOAL := ktss-extended

INIT = python3 init.py

SAVE = /opt/UPV/TFG/src/example/
VCF = -vcf $(SAVE)datosR1.vcf
FASTA = -fasta $(SAVE)hg19.fa.gz

LENGTH = 20
LENGTH_SEQUENCE = -p_p $(LENGTH) -p_s $(LENGTH)

KTSS_EXTENDED = -o pm -p e -m ktss
K = 4
KTSS_PARAMETERS = -k $(K) -ktss_nas False

RATIO = 0.90

test:
	poetry run pytest
black:
	poetry run black .
parser-extended:
	$(INIT) -p e $(VCF) $(FASTA) -s $(SAVE)
ktss-extended:
	$(INIT) $(KTSS_EXTENDED) $(VCF) $(FASTA) -s $(SAVE) $(LENGTH_SEQUENCE) -r $(RATIO) $(KTSS_PARAMETERS)
