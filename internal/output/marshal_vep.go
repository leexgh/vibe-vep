package output

import (
	"encoding/json"
	"strconv"

	"github.com/inodb/vibe-vep/internal/annotate"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// MarshalVEPAnnotation builds a VEPVariantAnnotation from a variant and its annotations,
// then marshals it to JSON. This is the standalone equivalent of JSONLWriter.marshalVEP().
func MarshalVEPAnnotation(input string, v *vcf.Variant, anns []*annotate.Annotation, assembly string) ([]byte, error) {
	ref, alt := alleleStrings(v)

	start, end := vepCoordinates(v)

	result := VEPVariantAnnotation{
		Input:         input,
		ID:            annotate.FormatVariantID(v.Chrom, v.Pos, v.Ref, v.Alt),
		SeqRegionName: v.Chrom,
		Start:         start,
		End:           end,
		AlleleString:  ref + "/" + alt,
		Strand:        1,
		AssemblyName:  assembly,
	}

	// Find most severe consequence.
	bestImpact := -1
	for _, ann := range anns {
		impact := annotate.ImpactRank(ann.Impact)
		if impact > bestImpact {
			bestImpact = impact
			result.MostSevereConsequence = firstConsequence(ann.Consequence)
		}
	}

	for _, ann := range anns {
		pStart, pEnd := proteinRange(ann)
		tc := VEPTranscriptConsequence{
			TranscriptID:        stripVersion(ann.TranscriptID),
			GeneID:              ann.GeneID,
			GeneSymbol:          ann.GeneName,
			GeneSymbolSource:    "HGNC",
			Biotype:             ann.Biotype,
			ConsequenceTerms:    splitConsequence(ann.Consequence),
			Impact:              ann.Impact,
			VariantAllele:       ann.Allele,
			AminoAcids:          formatAminoAcidsVEP(ann.AminoAcidChange),
			Codons:              ann.CodonChange,
			ProteinStart:        pStart,
			ProteinEnd:          pEnd,
			CDSStart:            ann.CDSPosition,
			CDSEnd:              ann.CDSPosition,
			CDNAStart:           ann.CDNAPosition,
			CDNAEnd:             ann.CDNAPosition,
			HGVSc:               prependTranscriptID(ann.TranscriptID, ann.HGVSc),
			HGVSp:               ann.HGVSp,
			Exon:                ann.ExonNumber,
			Intron:              ann.IntronNumber,
			RefSeqTranscriptIDs: ann.RefSeqIDs,
			HGVSOffset:          ann.HGVSOffset,
		}

		// SIFT/PolyPhen from annotation source extras.
		if s := ann.GetExtraKey("sift.score"); s != "" {
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				tc.SIFTScore = &f
			}
		}
		tc.SIFTPrediction = ann.GetExtraKey("sift.prediction")
		if s := ann.GetExtraKey("polyphen.score"); s != "" {
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				tc.PolyPhenScore = &f
			}
		}
		tc.PolyPhenPrediction = ann.GetExtraKey("polyphen.prediction")

		result.TranscriptConsequences = append(result.TranscriptConsequences, tc)
	}

	return json.Marshal(result)
}

// proteinRange returns the protein_start/protein_end VEP reports for an
// annotation. These are the UNSHIFTED positions of the change, which differ
// from the 3'-shifted position the HGVSp is written at: VEP reports PBRM1 as
// protein_start 1169 alongside its own p.I1170Sfs*23. Falls back to the HGVS
// position when the unshifted one was not computed (non-coding paths).
func proteinRange(ann *annotate.Annotation) (start, end int64) {
	start = ann.ProteinStart
	if start == 0 {
		start = ann.ProteinPosition
	}
	end = ann.ProteinEnd
	if end == 0 {
		end = start
	}
	return start, end
}
