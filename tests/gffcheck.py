def check_gff(path, chroms: dict) -> list:
    """Assert GFF3 validity: 9 columns, coordinates in range, unique IDs, Parents defined earlier and
    containing their children. Returns [(chrom, start, end, type, attrs)]."""
    ids, feats = {}, []
    with open(path) as f:
        assert f.readline().startswith("##gff-version 3")
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            assert len(x) == 9, line
            c, s, e = x[0], int(x[3]), int(x[4])
            assert c in chroms and 1 <= s <= e <= chroms[c], line
            a = dict(kv.split("=", 1) for kv in x[8].split(";"))
            assert a["ID"] not in ids, a["ID"]
            if "Parent" in a:
                pc, ps, pe = ids[a["Parent"]]
                assert pc == c and ps <= s and e <= pe, line
            ids[a["ID"]] = (c, s, e)
            feats.append((c, s, e, x[2], a))
    return feats
