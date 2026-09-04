package annotate

import "testing"

// Expected values are real Ensembl VEP output, taken from the genome-nexus
// test fixtures under service/src/test/resources/variant/.
func TestFormatCodonChangeCDS(t *testing.T) {
	tests := []struct {
		name             string
		cds              string
		delStart, delEnd int64
		inserted         string
		want             string
	}{
		// c.5946del -> agT/ag : deleted base is the last of its codon
		{"del at codon end", "aaaagt", 6, 6, "", "agT/ag"},
		// c.1856del -> cCt/ct : deleted base mid-codon
		{"del mid codon", "aaacct", 5, 5, "", "cCt/ct"},
		// c.1936del -> Agg/gg : deleted base first of codon
		{"del at codon start", "aaaagg", 4, 4, "", "Agg/gg"},
		// c.4049_4051del -> gAGGcc/gcc : three bases spanning two codons
		{"del spanning two codons", "aaagaggcc", 5, 7, "", "gAGGcc/gcc"},
		// c.29_30insCAG -> cag/caCAGg : insertion inside a codon
		{"insertion inside codon", "aaaaaacag", 9, 8, "CAG", "cag/caCAGg"},
		// c.627delinsAT -> atG/atAT
		{"delins one base to two", "aaaatg", 6, 6, "AT", "atG/atAT"},
		// c.1132_1133delinsT -> AGt/Tt
		{"delins two bases to one", "aaaagt", 4, 5, "T", "AGt/Tt"},
		// codon-boundary insertion disturbs no codon
		{"codon-boundary insertion", "aaaaaa", 7, 6, "CAACTT", "-/CAACTT"},
	}
	for _, tc := range tests {
		got := formatCodonChangeCDS(tc.cds, tc.delStart, tc.delEnd, tc.inserted)
		if got != tc.want {
			t.Errorf("%s: got %q, want %q", tc.name, got, tc.want)
		}
	}
}

func TestFormatCodonChangeMNV(t *testing.T) {
	tests := []struct {
		name       string
		cds        string
		start, end int64
		alt        string
		want       string
	}{
		{"DNP within one codon", "aaaggt", 4, 5, "AA", "GGt/AAt"},
		{"DNP spanning two codons", "aaaggtcc", 6, 7, "AA", "ggTCc/ggAAc"},
		{"single base", "aaaggt", 5, 5, "A", "gGt/gAt"},
	}
	for _, tc := range tests {
		got := formatCodonChangeMNV(tc.cds, tc.start, tc.end, tc.alt)
		if got != tc.want {
			t.Errorf("%s: got %q, want %q", tc.name, got, tc.want)
		}
	}
}
