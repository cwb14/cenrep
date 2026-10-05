import pytest

from cenrep.genomes import count_in, discover_genomes, load_fai, load_genes
from cenrep.util import CenrepError


def test_load_fai_missing(tmp_path):
    fa = tmp_path / "x.fa"
    fa.write_text(">c\nACGT\n")
    with pytest.raises(CenrepError, match="missing FASTA index"):
        load_fai(fa)


def _fa(d, name, chroms):
    (d / f"{name}.fa").write_text("".join(f">{c}\nACGT\n" for c in chroms))
    (d / f"{name}.fa.fai").write_text("".join(f"{c}\t4\t{i}\t4\t5\n" for i, c in enumerate(chroms)))


def test_discover_requires_a_genome(tmp_path):
    with pytest.raises(CenrepError, match="no \\*.fa genomes"):
        discover_genomes(tmp_path)


def test_discover_loads_any_number_of_genomes(tmp_path):
    _fa(tmp_path, "one", ["c1"])
    assert [g.name for g in discover_genomes(tmp_path)] == ["one"]
    for n in "bcde":
        _fa(tmp_path, n, ["c1", "c2"])
    gs = discover_genomes(tmp_path)
    assert [g.name for g in gs] == ["b", "c", "d", "e", "one"] and gs[0].chroms == {"c1": 4, "c2": 4}


def test_load_genes_reads_gene_gff(tmp_path):
    """genes/<genome>.gff3[.gz]: one midpoint per coding gene locus (isoforms counted once), native seqids,
    only sequences of the genome."""
    (tmp_path / "genes").mkdir()
    (tmp_path / "genes" / "Sp_A.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tx\tgene\t101\t200\t.\t+\t.\tID=g1;locus_type=coding\n"
        "chr1\tx\tmRNA\t101\t200\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "chr1\tx\tmRNA\t101\t180\t.\t+\t.\tID=g1.t2;Parent=g1\n"
        "chr1\tx\tgene\t1001\t1010\t.\t-\t.\tID=g2\n"
        "chr1\tx\tgene\t2001\t2010\t.\t-\t.\tID=g3;locus_type=no_orf\n"
        "scaf9\tx\tgene\t1\t10\t.\t+\t.\tID=g4;locus_type=coding\n")
    g = load_genes(tmp_path, "Sp_A", {"chr1": 5000, "chr2": 5000})
    assert set(g) == {"chr1"} and g["chr1"].tolist() == [150, 1005]
    assert set(load_genes(tmp_path, "Sp_A")) == {"chr1", "scaf9"}          # no chroms given: no filter
    assert load_genes(tmp_path, "Other", {"chr1": 5000}) is None


def test_load_genes_keeps_protein_coding_loci_of_typed_gffs(tmp_path):
    """NCBI/Ensembl `gene` lines include tRNA, ncRNA and pseudogene loci; with a biotype only protein_coding counts."""
    (tmp_path / "genes").mkdir()
    (tmp_path / "genes" / "Sp_A.gff3").write_text(
        "chr1\tRefSeq\tgene\t11\t20\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n"
        "chr1\tRefSeq\tgene\t101\t110\t.\t+\t.\tID=g2;gene_biotype=tRNA\n"
        "chr1\tEnsembl\tgene\t201\t210\t.\t+\t.\tID=g3;biotype=pseudogene\n"
        "chr1\tEnsembl\tgene\t301\t310\t.\t+\t.\tID=g4;biotype=protein_coding\n")
    assert load_genes(tmp_path, "Sp_A", {"chr1": 1000})["chr1"].tolist() == [15, 305]


def test_gene_gff_with_foreign_sequence_names_is_no_evidence(tmp_path):
    """Review I1: a genes/ file whose seqids match no sequence must not read as 'zero genes everywhere'."""
    (tmp_path / "genes").mkdir()
    (tmp_path / "genes" / "Sp_A.gff3").write_text(
        "Chr01\tx\tgene\t1\t10\t.\t+\t.\tID=g1;orig_locus_type=other\n")
    assert load_genes(tmp_path, "Sp_A", {"chr1": 5000}) is None
    assert load_genes(tmp_path, "Sp_A", {"Chr01": 5000})["Chr01"].tolist() == [5]   # orig_locus_type is not locus_type


def test_gene_gff_bad_coordinates_fail_clearly(tmp_path):
    (tmp_path / "genes").mkdir()
    (tmp_path / "genes" / "Sp_A.gff3").write_text("chr1\tx\tgene\tone\t10\t.\t+\t.\tID=g1\n")
    with pytest.raises(CenrepError, match="Sp_A.gff3:1"):
        load_genes(tmp_path, "Sp_A", {"chr1": 5000})


def test_load_genes_reads_gzipped_gene_gff(tmp_path):
    import gzip
    (tmp_path / "genes").mkdir()
    with gzip.open(tmp_path / "genes" / "Sp_A.gff3.gz", "wt") as f:
        f.write("chr1\tx\tgene\t11\t20\t.\t+\t.\tID=g1\n")
    assert load_genes(tmp_path, "Sp_A", {"chr1": 100})["chr1"].tolist() == [15]


def _bgzip_fa(d, name, chroms):
    """Bgzip-compressed FASTA with .fai + .gzi, as the cleaned polyploid genomes are stored."""
    import pysam
    plain = d / f"{name}.plain.fa"
    plain.write_text("".join(f">{c} orig=x\n{'ACGT' * 30}\n" for c in chroms))
    gz = d / f"{name}.fa.gz"
    pysam.tabix_compress(str(plain), str(gz), force=True)
    plain.unlink()
    pysam.faidx(str(gz))
    return gz


def test_discover_reads_bgzipped_genomes(tmp_path):
    import pysam
    gz = _bgzip_fa(tmp_path, "Sp.2x.progenitor-A.GCA_1.1", ["c1", "c2"])
    _fa(tmp_path, "Other", ["x1"])
    gs = discover_genomes(tmp_path)
    assert [g.name for g in gs] == ["Other", "Sp.2x.progenitor-A.GCA_1.1"]
    g = gs[1]
    assert g.fasta == gz and g.chroms == {"c1": 120, "c2": 120}
    with pysam.FastaFile(str(g.fasta)) as fa:
        assert fa.fetch("c2", 4, 12) == "ACGTACGT"


def test_discover_rejects_gz_without_bgzip_index(tmp_path):
    gz = _bgzip_fa(tmp_path, "Sp_A", ["c1"])
    (tmp_path / "Sp_A.fa.gz.gzi").unlink()
    with pytest.raises(CenrepError, match="bgzip"):
        discover_genomes(tmp_path)


def test_discover_rejects_one_genome_twice(tmp_path):
    _bgzip_fa(tmp_path, "Sp_A", ["c1"])
    _fa(tmp_path, "Sp_A", ["c1"])
    with pytest.raises(CenrepError, match="Sp_A"):
        discover_genomes(tmp_path)



def test_discover_reads_fasta_and_fna_names(tmp_path):
    for n in ("A.fasta", "B.fna", "C.fa"):
        (tmp_path / n).write_text(">c1\nACGT\n")
        (tmp_path / f"{n}.fai").write_text("c1\t4\t4\t4\t5\n")
    (tmp_path / "notes.txt").write_text("x")
    assert [g.name for g in discover_genomes(tmp_path)] == ["A", "B", "C"]
