from src.vcf.vcfReader import VcfMutationsReader

results = [
    ('', 'T', 'ATCAA'),
    ('T', 'A', 'TCAAT'),
    ('TA', 'T', 'CAATG'),
    ('TAT', 'C', 'AATGC'),
    ('TATC', 'A', 'ATGCC'),
    ('TATCA', 'A', 'TGCCT'),
    ('GAGGG', 'C', 'ATAAC'),
    ('AGGGC', 'A', 'TAACA'),
    ('GGGCA', 'T', 'AACAT'),
    ('GGCAT', 'A', 'ACATC'),
    ('AGGTA', 'T', 'CTAAT'),
    ('GGTAT', 'C', 'TAAT'),
    ('GTATC', 'T', 'AAT'),
    ('TATCT', 'A', 'AT'),
    ('ATCTA', 'A', 'T'),
    ('TCTAA', 'T', '')
]


def test_get_seq_by_chr_pos():
    reader = VcfMutationsReader(
        '/opt/UPV/TFG/src/tests/vcfTest.vcf',
        '/opt/UPV/TFG/src/tests/test.fa.gz'
    )

    for test_vcf_line in zip(reader.get_vcf(), results * 3):
        vcf_line = test_vcf_line[0]
        result = test_vcf_line[1]
        reference_seq_fasta = reader.get_seq_by_chr_pos(vcf_line.CHROM, vcf_line.POS, 5, 5)

        assert result[0] == reference_seq_fasta[0]
        assert vcf_line.REF == reference_seq_fasta[1]
        assert result[2] == reference_seq_fasta[2]
