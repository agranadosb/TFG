test:
	poetry run pytest
black:
	poetry run black .
parser-extended:
	python3 init.py -p e -vcf /opt/UPV/TFG/src/example/datosR1.vcf -fasta /opt/UPV/TFG/src/example/hg19.fa.gz -s /opt/UPV/TFG/src/example/
model-ktss:
	python3 init.py -o m -vcf /opt/UPV/TFG/src/example/datosR1.vcf -fasta /opt/UPV/TFG/src/example/hg19.fa.gz -s /opt/UPV/TFG/src/example/
parser-extended-model-ktss:
	python3 init.py -o pm -p e -m ktss -vcf /opt/UPV/TFG/src/example/datosR1.vcf -fasta /opt/UPV/TFG/src/example/hg19.fa.gz -s /opt/UPV/TFG/src/example/
