from src.vcf.vcfReader import VcfMutationsReader


def test_get_seq_by_chr_pos():
    reader = VcfMutationsReader(
        '/opt/UPV/TFG/src/tests/vcfTest.vcf',
        '/opt/UPV/TFG/src/tests/test.fa.gz'
    )
    for i in reader.getVcf():
        reference = i.REF
        reference_seq_fasta = reader.get_seq_by_chr_pos(i.CHROM, i.POS, 5, 5)

        print(reference)
        print(reference_seq_fasta)
        print()
        assert reference == reference_seq_fasta[1]
