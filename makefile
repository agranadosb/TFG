.DEFAULT_GOAL := ktss-extended

INIT = python3 init.py

SAVE = /opt/UPV/TFG/src/example/
VCF = -vcf $(SAVE)datosR1.vcf
FASTA = -fasta $(SAVE)hg19.fa.gz

LENGTH_SEQUENCE = 20
PS = -p_p $(LENGTH_SEQUENCE) -p_s $(LENGTH_SEQUENCE)

KTSS_EXTENDED = -o pm -p e -m ktss

RATIO = -r 0.90

test:
	poetry run pytest
black:
	poetry run black .
parser-extended:
	$(INIT) -p e $(VCF) $(FASTA) -s $(SAVE)
ktss-extended:
	$(INIT) $(KTSS_EXTENDED) $(VCF) $(FASTA) -s $(SAVE) $(PS) $(RATIO)
