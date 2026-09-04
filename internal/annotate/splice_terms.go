package annotate

import "github.com/inodb/vibe-vep/internal/cache"

// Additional SO terms VEP reports for positions inside an intron near a splice
// site. vibe-vep previously collapsed all of these into splice_region_variant
// or intron_variant.
const (
	ConsequenceSpliceDonor5thBase  = "splice_donor_5th_base_variant"
	ConsequenceSpliceDonorRegion   = "splice_donor_region_variant"
	ConsequencePolypyrimidineTract = "splice_polypyrimidine_tract_variant"
	ConsequenceProteinAltering     = "protein_altering_variant"
)

// consequenceRank orders terms the way Ensembl ranks consequence severity,
// which is the order VEP emits them in. Only the terms vibe-vep can produce
// need an entry; anything unknown sorts last but keeps a stable relative order.
var consequenceRank = map[string]int{
	ConsequenceSpliceAcceptor:        2,
	ConsequenceSpliceDonor:           3,
	ConsequenceStopGained:            4,
	ConsequenceFrameshiftVariant:     5,
	ConsequenceStopLost:              6,
	ConsequenceStartLost:             7,
	ConsequenceInframeInsertion:      9,
	ConsequenceInframeDeletion:       10,
	ConsequenceMissenseVariant:       11,
	ConsequenceProteinAltering:       12,
	ConsequenceSpliceDonor5thBase:    13,
	ConsequenceSpliceRegion:          14,
	ConsequenceSpliceDonorRegion:     15,
	ConsequencePolypyrimidineTract:   16,
	ConsequenceStartRetained:         18,
	ConsequenceStopRetained:          19,
	ConsequenceSynonymousVariant:     20,
	ConsequenceCodingSequenceVariant: 21,
	ConsequenceMatureMiRNA:           22,
	Consequence5PrimeUTR:             23,
	Consequence3PrimeUTR:             24,
	ConsequenceNonCodingExon:         25,
	ConsequenceIntronVariant:         26,
	ConsequenceUpstreamGene:          29,
	ConsequenceDownstreamGene:        30,
	ConsequenceIntergenicVariant:     31,
}

func rankOf(term string) int {
	if r, ok := consequenceRank[term]; ok {
		return r
	}
	return 100
}

// intronOffset returns the signed distance from pos to the nearest exon
// boundary, in transcript orientation: positive N means the Nth base into the
// intron on the donor side, negative N the Nth base before the next exon on the
// acceptor side. ok is false when pos is not intronic for this transcript.
func intronOffset(pos int64, t *cache.Transcript) (offset int64, ok bool) {
	idx := t.FindNearestExonIdx(pos)
	if idx < 0 {
		return 0, false
	}
	best := int64(0)
	found := false
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= len(t.Exons) {
			continue
		}
		exon := &t.Exons[i]
		if pos >= exon.Start && pos <= exon.End {
			return 0, false // exonic
		}
		var d int64
		if pos > exon.End {
			d = pos - exon.End // downstream of this exon
			if !t.IsForwardStrand() {
				d = -d
			}
		} else {
			d = exon.Start - pos // upstream of this exon
			if t.IsForwardStrand() {
				d = -d
			}
		}
		if !found || abs64(d) < abs64(best) {
			best, found = d, true
		}
	}
	return best, found
}

// spliceTermsForOffset returns the SO terms VEP reports for a single intronic
// position at the given signed offset.
//
// Derived from single-base substitutions in a VEP111 MSK-IMPACT MAF:
//
//	+1,+2   splice_donor_variant                                  (no intron_variant)
//	+3,+4,+6 splice_donor_region_variant, intron_variant
//	+5      splice_donor_5th_base_variant, intron_variant         (not donor_region)
//	+7,+8   splice_region_variant, intron_variant
//	-1,-2   splice_acceptor_variant                               (no intron_variant)
//	-3..-8  splice_region_variant, splice_polypyrimidine_tract_variant, intron_variant
//	-9..-17 splice_polypyrimidine_tract_variant, intron_variant
//	beyond  intron_variant
func spliceTermsForOffset(offset int64) []string {
	switch {
	case offset == 1 || offset == 2:
		return []string{ConsequenceSpliceDonor}
	case offset == 5:
		return []string{ConsequenceSpliceDonor5thBase, ConsequenceIntronVariant}
	case offset >= 3 && offset <= 6:
		return []string{ConsequenceSpliceDonorRegion, ConsequenceIntronVariant}
	case offset >= 7 && offset <= 8:
		return []string{ConsequenceSpliceRegion, ConsequenceIntronVariant}
	case offset == -1 || offset == -2:
		return []string{ConsequenceSpliceAcceptor}
	case offset <= -3 && offset >= -8:
		return []string{ConsequenceSpliceRegion, ConsequencePolypyrimidineTract, ConsequenceIntronVariant}
	case offset <= -9 && offset >= -17:
		return []string{ConsequencePolypyrimidineTract, ConsequenceIntronVariant}
	default:
		return []string{ConsequenceIntronVariant}
	}
}

// joinConsequenceTerms renders a term set in Ensembl severity order, dropping
// duplicates. splice_region_variant is suppressed when a splice donor or
// acceptor term is present, which is what VEP does.
func joinConsequenceTerms(terms []string) string {
	seen := make(map[string]bool, len(terms))
	uniq := make([]string, 0, len(terms))
	hasSite := false
	has5thBase := false
	for _, t := range terms {
		if t == "" || seen[t] {
			continue
		}
		switch t {
		case ConsequenceSpliceDonor, ConsequenceSpliceAcceptor:
			hasSite = true
		case ConsequenceSpliceDonor5thBase:
			has5thBase = true
		}
		seen[t] = true
		uniq = append(uniq, t)
	}
	out := make([]string, 0, len(uniq))
	for _, t := range uniq {
		// splice_region_variant is not reported alongside a donor/acceptor
		// site, and splice_donor_5th_base_variant is the more specific term
		// where it overlaps splice_donor_region_variant. Both suppressions are
		// what VEP emits for spans covering these positions.
		if hasSite && t == ConsequenceSpliceRegion {
			continue
		}
		if has5thBase && t == ConsequenceSpliceDonorRegion {
			continue
		}
		out = append(out, t)
	}
	// insertion sort by rank; the lists are tiny
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && rankOf(out[j]) < rankOf(out[j-1]); j-- {
			out[j], out[j-1] = out[j-1], out[j]
		}
	}
	res := ""
	for i, t := range out {
		if i > 0 {
			res += ","
		}
		res += t
	}
	return res
}
