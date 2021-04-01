from src.model.parserVcf import ParserVcf


def test_generate_simplified_sequences():
    parser = ParserVcf(
        '/opt/UPV/TFG/src/tests/vcfTest.vcf',
        '/opt/UPV/TFG/src/tests/test.fa.gz'
    )

    # TODO: Hay que testear get_simplified_sequence

    parser.generate_simplified_sequences('/opt/UPV/TFG/src/tests', write_chromosme=True)
