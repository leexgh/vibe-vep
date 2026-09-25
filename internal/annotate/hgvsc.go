package annotate

import (
	"strconv"

	"github.com/inodb/vibe-vep/internal/cache"
	"github.com/inodb/vibe-vep/internal/vcf"
)

// FormatHGVSc formats the HGVS coding DNA notation for a variant on a transcript.
// Returns empty string for non-coding transcripts or upstream/downstream variants.
//
// The function is optimized for minimal allocations on CDS-based paths (SNV, deletion,
// insertion) by using stack-allocated buffers and lazy computation of reverse complement.
func FormatHGVSc(v *vcf.Variant, t *cache.Transcript, result *ConsequenceResult) string {
	if t == nil || v == nil || result == nil {
		return ""
	}

	// Skip upstream/downstream/intergenic
	switch result.Consequence {
	case ConsequenceUpstreamGene, ConsequenceDownstreamGene, ConsequenceIntergenicVariant:
		return ""
	}

	// Determine prefix: c. for coding, n. for non-coding
	prefix := "c."
	if !t.IsProteinCoding() {
		prefix = "n."
	}

	// SNV and MNV (same length ref/alt): compute RC and format.
	if !v.IsIndel() {
		// No change (ref==alt): HGVS uses "=" notation
		if v.Ref == v.Alt {
			startPosStr := genomicToHGVScPos(v.Pos, t)
			if startPosStr == "" {
				return ""
			}
			return prefix + startPosStr + "="
		}
		ref := v.Ref
		alt := v.Alt
		if t.IsReverseStrand() {
			ref = ReverseComplement(ref)
			alt = ReverseComplement(alt)
		}
		startPosStr := genomicToHGVScPos(v.Pos, t)
		if startPosStr == "" {
			return ""
		}
		// SNV: single nucleotide substitution (c.123A>G)
		if len(v.Ref) == 1 {
			return prefix + startPosStr + ref + ">" + alt
		}
		// MNV: multi-nucleotide changes must use delins per HGVS
		// (substitution is defined as exactly one nucleotide replaced by one)
		endPosStr := genomicToHGVScPos(v.Pos+int64(len(v.Ref))-1, t)
		if t.IsReverseStrand() {
			startPosStr, endPosStr = endPosStr, startPosStr
		}
		if endPosStr == "" {
			return ""
		}
		return prefix + startPosStr + "_" + endPosStr + "delins" + alt
	}

	// Indel handling — defer RC computation to branches that need it.
	refLen := len(v.Ref)
	altLen := len(v.Alt)

	if altLen > refLen {
		return formatHGVScInsertion(v, t, prefix, refLen, altLen, result)
	}

	return formatHGVScDeletion(v, t, prefix, refLen, altLen, result)
}

// formatHGVScInsertion handles the insertion path of FormatHGVSc.
// It computes the inserted sequence on the coding strand using a stack buffer
// to avoid allocations from ReverseComplement.
func formatHGVScInsertion(v *vcf.Variant, t *cache.Transcript, prefix string, refLen, altLen int, result *ConsequenceResult) string {
	insLen := altLen - refLen

	// Compute inserted sequence on coding strand into a stack buffer.
	// For reverse strand, we write the reverse complement directly to avoid
	// allocating a string from ReverseComplement().
	var seqStack [64]byte
	var seq []byte
	if insLen <= len(seqStack) {
		seq = seqStack[:insLen]
	} else {
		seq = make([]byte, insLen)
	}

	if t.IsReverseStrand() {
		// RC(v.Alt)[:altLen-refLen] == RC(v.Alt[refLen:])
		genomicInserted := v.Alt[refLen:]
		for i := 0; i < insLen; i++ {
			seq[i] = Complement(genomicInserted[insLen-1-i])
		}
	} else {
		copy(seq, v.Alt[refLen:])
	}

	// 3' shift insertion in CDS space, then check dup at shifted position
	cdsPos := GenomicToCDS(v.Pos, t)
	if cdsPos > 0 && len(t.CDSSequence) > 0 {
		cdsIdx := int(cdsPos - 1) // 0-based anchor index
		if t.IsReverseStrand() {
			cdsIdx-- // RC makes shared base a suffix; insertion is before CDS position
		}
		if cdsIdx < 0 {
			cdsIdx = 0
		}

		// Shift insertion in place (modifies seq), stops at exon boundary
		maxIdx := cdsExonEndIdx(cdsIdx, t)
		shiftedIdx := shiftInsertionBuf(seq, cdsIdx, t.CDSSequence, maxIdx)
		result.HGVSOffset = shiftedIdx - cdsIdx
		// Codons are reported at the UNSHIFTED position. VEP writes the HGVS at
		// the 3'-shifted position but expands codons around where the change
		// actually sits, so subtract the shift back out.
		unshiftedAnchor := int64(shiftedIdx+1) - int64(result.HGVSOffset)
		result.CodonChange = formatCodonChangeCDS(t.CDSSequence,
			unshiftedAnchor+1, unshiftedAnchor, string(seq))
		seqLen := len(seq)

		// HGVS defines insertion and duplication as mutually exclusive: an
		// insertion is a change "where the insertion is not a copy of a sequence
		// immediately 5'", and a duplication is a copy inserted "directly 3' of
		// the original copy of that sequence". So whenever the inserted bases are
		// a tandem copy of what precedes them, "ins" is not a valid description.
		//
		// This deliberately diverges from Ensembl VEP, which only writes "dup"
		// when the 3' shift moved the variant. Across a VEP111 MSK-IMPACT MAF
		// 6,769 of VEP's 8,124 "ins" descriptions are tandem copies by the
		// reference sequence -- e.g. c.1663_1664insC inserted into a run of seven
		// C's, which HGVS requires be written c.1663dup. VEP also contradicts
		// itself on 541 rows, writing "ins" at the CDS level and "dup" at the
		// protein level for the same variant (ASXL2 c.353_354insCAG / p.S117dup).
		// Matching VEP here would cost nothing in accuracy but would emit
		// descriptions that fail to match ClinVar and LOVD, which normalize to
		// "dup". Concordance with vep111 on HGVSc is ~1.4pp lower as a result.

		// Check dup: inserted bases match preceding bases at shifted position.
		// Go optimizes string([]byte) == string comparisons to avoid allocation.
		dupStart := shiftedIdx - seqLen + 1
		if dupStart >= 0 && shiftedIdx+1 <= len(t.CDSSequence) &&
			t.CDSSequence[dupStart:shiftedIdx+1] == string(seq) {
			return cdsPosRangeStr(dupStart+1, shiftedIdx+1, "dup")
		}

		// Check dup: inserted bases match following bases at shifted position
		afterStart := shiftedIdx + 1
		afterEnd := afterStart + seqLen
		if afterStart >= 0 && afterEnd <= len(t.CDSSequence) &&
			t.CDSSequence[afterStart:afterEnd] == string(seq) {
			return cdsPosRangeStr(afterStart+1, afterEnd, "dup")
		}

		// Plain insertion at shifted CDS position
		if shiftedIdx+2 <= len(t.CDSSequence) {
			var buf [128]byte
			n := copy(buf[:], "c.")
			n += putInt64(buf[n:], int64(shiftedIdx+1))
			buf[n] = '_'
			n++
			n += putInt64(buf[n:], int64(shiftedIdx+2))
			n += copy(buf[n:], "ins")
			n += copy(buf[n:], seq)
			return string(buf[:n])
		}
		// Insertion at end of CDS: second flanking position is *1 (3'UTR)
		if shiftedIdx+1 == len(t.CDSSequence) {
			return prefix + strconv.FormatInt(int64(shiftedIdx+1), 10) + "_*1ins" + string(seq)
		}
	}

	// Fall back to genomic-based positions for non-CDS insertions.
	// These paths use string concat (less common, not on the hot path).
	insertedSeq := string(seq)

	dup := checkDuplication(v, t, insertedSeq)
	if dup.isDup {
		if dup.cdsStart == dup.cdsEnd {
			return "c." + strconv.FormatInt(dup.cdsStart, 10) + "dup"
		}
		return "c." + strconv.FormatInt(dup.cdsStart, 10) + "_" + strconv.FormatInt(dup.cdsEnd, 10) + "dup"
	}

	// Check for duplication at splice junctions where one flanking position
	// is intronic (CDS=0) and the other is exonic. The CDS-based dup check
	// above bails when the anchor is intronic, but the inserted base may
	// still duplicate the exon boundary base.
	if sjDup := checkSpliceJunctionDup(v, t, insertedSeq); sjDup.isDup {
		if sjDup.cdsStart == sjDup.cdsEnd {
			return "c." + strconv.FormatInt(sjDup.cdsStart, 10) + "dup"
		}
		return "c." + strconv.FormatInt(sjDup.cdsStart, 10) + "_" + strconv.FormatInt(sjDup.cdsEnd, 10) + "dup"
	}

	insAfterPos := v.Pos
	insBeforePos := v.Pos + 1
	if t.IsReverseStrand() {
		insAfterPos = v.Pos + 1
		insBeforePos = v.Pos
	}
	afterStr := genomicToHGVScPos(insAfterPos, t)
	beforeStr := genomicToHGVScPos(insBeforePos, t)
	return prefix + afterStr + "_" + beforeStr + "ins" + insertedSeq
}

// formatHGVScDeletion handles the deletion (and delins) path of FormatHGVSc.
// It avoids calling ReverseComplement for the ref/alt at the top since the
// deletion CDS path only needs CDS positions and possibly RC of the extra alt bases.
func formatHGVScDeletion(v *vcf.Variant, t *cache.Transcript, prefix string, refLen, altLen int, result *ConsequenceResult) string {
	// Compute actual shared prefix length on genomic strand
	sharedLen := 0
	for sharedLen < refLen && sharedLen < altLen && v.Ref[sharedLen] == v.Alt[sharedLen] {
		sharedLen++
	}
	if sharedLen == 0 {
		sharedLen = 1
	}
	if sharedLen > altLen {
		sharedLen = altLen
	}

	delStartGenomic := v.Pos + int64(sharedLen)
	delEndGenomic := v.Pos + int64(refLen) - 1
	extraAlt := ""
	if sharedLen < altLen {
		extraAlt = v.Alt[sharedLen:]
	}

	// Clip shared suffix between remaining ref and extra alt.
	if len(extraAlt) > 0 {
		remainRef := v.Ref[sharedLen:]
		sharedSuffix := 0
		for sharedSuffix < len(remainRef) && sharedSuffix < len(extraAlt) &&
			remainRef[len(remainRef)-1-sharedSuffix] == extraAlt[len(extraAlt)-1-sharedSuffix] {
			sharedSuffix++
		}
		if sharedSuffix > 0 {
			delEndGenomic -= int64(sharedSuffix)
			extraAlt = extraAlt[:len(extraAlt)-sharedSuffix]
		}
	}

	if t.IsReverseStrand() {
		delStartGenomic, delEndGenomic = delEndGenomic, delStartGenomic
	}

	// Try CDS-based positioning (hot path — optimized with stack buffers)
	delStartCDS := GenomicToCDS(delStartGenomic, t)
	delEndCDS := GenomicToCDS(delEndGenomic, t)
	if delStartCDS > 0 && delEndCDS > 0 && len(t.CDSSequence) > 0 {
		if len(extraAlt) > 0 {
			// Delins: no 3' shift per HGVS convention.
			codingExtra := extraAlt
			if t.IsReverseStrand() {
				codingExtra = ReverseComplement(extraAlt)
			}
			result.CodonChange = formatCodonChangeCDS(t.CDSSequence,
				delStartCDS, delEndCDS, codingExtra)
			// Write RC of extra alt directly into the output buffer for reverse strand.
			var buf [128]byte
			n := copy(buf[:], "c.")
			n += putInt64(buf[n:], delStartCDS)
			if delStartCDS != delEndCDS {
				buf[n] = '_'
				n++
				n += putInt64(buf[n:], delEndCDS)
			}
			n += copy(buf[n:], "delins")
			if t.IsReverseStrand() {
				// Write reverse complement of extraAlt directly into buffer
				eLen := len(extraAlt)
				for i := 0; i < eLen; i++ {
					buf[n+i] = Complement(extraAlt[eLen-1-i])
				}
				n += eLen
			} else {
				n += copy(buf[n:], extraAlt)
			}
			return string(buf[:n])
		}
		// A deletion whose 3' end reaches the exon edge may keep shifting into
		// the intron (see deletionShiftIntoIntron).
		if sc, ee, _, k, ok := deletionShiftIntoIntron(v, t); ok {
			result.HGVSOffset = int(sc-delStartCDS) + 0
			if result.HGVSOffset < 0 {
				result.HGVSOffset = 0
			}
			return formatDeletionAcrossDonor(sc, ee, k)
		}

		// Pure deletion: apply 3' shift (stops at exon boundary)
		maxIdx := cdsExonEndIdx(int(delEndCDS-1), t)
		sStart, sEnd := shiftDeletionThreePrime(int(delStartCDS-1), int(delEndCDS-1), t.CDSSequence, maxIdx)
		result.HGVSOffset = sStart - int(delStartCDS-1)
		// Codons are expanded around the unshifted deleted range (see above).
		result.CodonChange = formatCodonChangeCDS(t.CDSSequence,
			int64(sStart+1-result.HGVSOffset), int64(sEnd+1-result.HGVSOffset), "")
		return cdsPosRangeStr(sStart+1, sEnd+1, "del")
	}

	// Before falling back, a purely intronic pure deletion may still 3'-shift
	// across the junction into the exon (see shiftedIntronicDelCDS).
	if len(extraAlt) == 0 {
		if sStart, sEnd, k := shiftedIntronicDelCDS(v, t); sStart > 0 {
			result.HGVSOffset = k
			return cdsPosRangeStr(int(sStart), int(sEnd), "del")
		}
	}

	// Fall back to genomic-based positions for intronic/UTR/non-coding deletions.
	// Convert extra alt bases to coding strand (only needed for fallback path).
	codingExtraAlt := extraAlt
	if t.IsReverseStrand() && len(extraAlt) > 0 {
		codingExtraAlt = ReverseComplement(extraAlt)
	}

	// A deletion running from an exon into the following intron can still be
	// 3'-shifted, using the stored intron flank (see straddlingDeletionShift).
	if len(extraAlt) == 0 {
		if k, _ := straddlingDeletionShift(v, t); k > 0 {
			d := int64(k)
			if t.IsReverseStrand() {
				d = -d
			}
			delStartGenomic += d
			delEndGenomic += d
			result.HGVSOffset = k
		}
	}

	delStartStr := genomicToHGVScPos(delStartGenomic, t)
	if len(codingExtraAlt) > 0 {
		if delStartGenomic == delEndGenomic {
			return prefix + delStartStr + "delins" + codingExtraAlt
		}
		delEndStr := genomicToHGVScPos(delEndGenomic, t)
		return prefix + delStartStr + "_" + delEndStr + "delins" + codingExtraAlt
	}
	if delStartGenomic == delEndGenomic {
		return prefix + delStartStr + "del"
	}
	delEndStr := genomicToHGVScPos(delEndGenomic, t)
	return prefix + delStartStr + "_" + delEndStr + "del"
}

// deletionShiftIntoIntron handles a pure deletion whose 3' end reaches the end
// of an exon and can keep shifting 3' into the following intron.
//
// HGVS places a deletion as far 3' as possible and does not stop at a splice
// boundary. CYLD c.2099del sits on the last base of an exon whose intron opens
// with the same base, so VEP reports c.2099+1del: the deletion has moved onto
// the essential donor base. vibe-vep stopped at the exon end and reported the
// unshifted c.2099del.
//
// Returns the post-shift exonic start (1-based CDS), the exon's last CDS
// position, the deletion length, and how many bases it moved into the intron.
// ok is false when the shape does not apply or no flank is stored.
func deletionShiftIntoIntron(v *vcf.Variant, t *cache.Transcript) (startCDS, exonEndCDS int64, delLen, k int, ok bool) {
	refLen := len(v.Ref)
	if refLen == 0 || len(v.Alt) != 0 || len(t.CDSSequence) == 0 || !t.IsProteinCoding() {
		return 0, 0, 0, 0, false
	}
	lo, hi := v.Pos, v.Pos+int64(refLen)-1
	codingStart, codingEnd := lo, hi
	if t.IsReverseStrand() {
		codingStart, codingEnd = hi, lo
	}
	sCDS := GenomicToCDS(codingStart, t)
	eCDS := GenomicToCDS(codingEnd, t)
	if sCDS < 1 || eCDS < 1 {
		return 0, 0, 0, 0, false // not fully exonic
	}
	exon := t.FindExon(codingEnd)
	if exon == nil {
		return 0, 0, 0, 0, false
	}
	exonCodingEnd := exon.End
	if t.IsReverseStrand() {
		exonCodingEnd = exon.Start
	}
	eEnd := GenomicToCDS(exonCodingEnd, t)
	if eEnd < 1 {
		return 0, 0, 0, 0, false
	}

	// Shift within the exon first, exactly as the ordinary path does.
	sStart, sEnd := shiftDeletionThreePrime(int(sCDS-1), int(eCDS-1), t.CDSSequence, int(eEnd-1))
	if int64(sEnd) != eEnd-1 {
		return 0, 0, 0, 0, false // 3' end never reaches the exon edge; nothing to cross
	}
	s := int64(sStart) + 1
	L := sEnd - sStart + 1

	// Keep going into the intron. Position p is CDS while p <= eEnd, and the
	// (p-eEnd)th intron base beyond that, so one accessor covers both sides.
	flank := exon.IntronAfter
	cds := t.CDSSequence
	baseAt := func(p int64) (byte, bool) {
		if p <= eEnd {
			if p < 1 || int(p) > len(cds) {
				return 0, false
			}
			return cds[p-1], true
		}
		i := int(p - eEnd - 1)
		if i >= len(flank) {
			return 0, false
		}
		return flank[i], true
	}
	n := 0
	for {
		leaving, okL := baseAt(s + int64(n))
		entering, okE := baseAt(eEnd + int64(n) + 1)
		if !okL || !okE || leaving != entering {
			break
		}
		n++
	}
	if n == 0 {
		return 0, 0, 0, 0, false
	}
	return s + int64(n), eEnd, L, n, true
}

// formatDeletionAcrossDonor renders a deletion that has shifted k bases past the
// end of an exon, keeping whatever exonic part remains.
func formatDeletionAcrossDonor(startCDS, exonEndCDS int64, k int) string {
	if startCDS <= exonEndCDS {
		// Part of the deletion is still exonic: c.START_EXONEND+k del
		return "c." + strconv.FormatInt(startCDS, 10) + "_" +
			strconv.FormatInt(exonEndCDS, 10) + "+" + strconv.Itoa(k) + "del"
	}
	// Entirely intronic now: c.EXONEND+a_EXONEND+k del
	a := int(startCDS - exonEndCDS)
	base := "c." + strconv.FormatInt(exonEndCDS, 10) + "+"
	if a == k {
		return base + strconv.Itoa(k) + "del"
	}
	return base + strconv.Itoa(a) + "_" + base + strconv.Itoa(k) + "del"
}

// deletionHitsDonor reports whether the shifted deletion covers the first two
// intron bases, the essential donor dinucleotide.
func deletionHitsDonor(delLen, k int) bool {
	first := k - delLen + 1
	if first < 1 {
		first = 1
	}
	return first <= 2 && k >= 1
}

// straddlingDeletionShift computes how far a pure deletion running from an exon
// into the following intron may be shifted 3'.
//
// HGVS places a deletion as far 3' as possible. Each step is allowed when the
// base leaving the front of the deleted span equals the base entering at the
// back; here the front is exonic and the back is intronic, so the test needs
// intron sequence. That is what Exon.IntronAfter carries.
//
//	MET  vibe c.3015_3028+6del -> VEP c.3016_3028+7del  (offset 1)
//	B2M  vibe c.59_67+237del   -> VEP c.61_67+239del    (offset 2)
//
// Returns the shift distance and the CDS position the deletion starts at after
// shifting. Returns 0 when the shape does not apply or no flank is stored, in
// which case the caller keeps the unshifted description.
func straddlingDeletionShift(v *vcf.Variant, t *cache.Transcript) (shift int, newStartCDS int64) {
	refLen := len(v.Ref)
	if refLen == 0 || len(v.Alt) != 0 || len(t.CDSSequence) == 0 || !t.IsProteinCoding() {
		return 0, 0
	}
	lo, hi := v.Pos, v.Pos+int64(refLen)-1

	// Orient the span: codingStart is the 5' end of the deletion on the coding
	// strand, codingEnd the 3' end.
	codingStart, codingEnd := lo, hi
	if t.IsReverseStrand() {
		codingStart, codingEnd = hi, lo
	}
	startCDS := GenomicToCDS(codingStart, t)
	if startCDS < 1 || GenomicToCDS(codingEnd, t) > 0 {
		return 0, 0 // not a straddle: either the 5' end is intronic or the 3' end is still coding
	}
	exon := t.FindExon(codingStart)
	if exon == nil {
		return 0, 0
	}

	// How many intronic bases the deletion already covers.
	exonCodingEnd := exon.End
	if t.IsReverseStrand() {
		exonCodingEnd = exon.Start
	}
	tail := codingEnd - exonCodingEnd
	if t.IsReverseStrand() {
		tail = exonCodingEnd - codingEnd
	}
	if tail < 1 {
		return 0, 0
	}

	flank := exon.IntronAfter
	cds := t.CDSSequence
	k := 0
	for {
		fi := int(tail) + k // 0-based index of the base just past the deletion
		ci := int(startCDS) - 1 + k
		if fi >= len(flank) || ci >= len(cds) || flank[fi] != cds[ci] {
			break
		}
		k++
	}
	if k == 0 {
		return 0, 0
	}
	return k, startCDS + int64(k)
}

// shiftedIntronicDelCDS handles a pure deletion lying entirely within an intron
// immediately 5' (on the coding strand) of a coding exon, where the deleted
// bases keep repeating into that exon. HGVS assigns the most 3' position
// possible, so such a deletion slides across the splice junction: POLD1 deletes
// a G from the last intron base of an exon that opens GGGGGG, and VEP reports
// c.2959del with hgvs_offset 6, not c.2954-1del.
//
// No intron sequence is needed, which matters because vibe-vep has none. The
// deleted bases come from the variant itself and the run continues in
// CDSSequence, so sliding k bases is allowed while the exon base at j matches
// deleted[j%L]. Once k reaches L the deletion lies wholly inside the exon;
// below that it still straddles the junction and HGVS keeps the intronic form.
//
// Returns the 1-based CDS range it lands on and the shift distance, or zeros.
func shiftedIntronicDelCDS(v *vcf.Variant, t *cache.Transcript) (cdsStart, cdsEnd int64, shift int) {
	refLen := len(v.Ref)
	if refLen == 0 || len(v.Alt) != 0 || len(t.CDSSequence) == 0 || !t.IsProteinCoding() {
		return 0, 0, 0
	}
	delStartGenomic := v.Pos
	delEndGenomic := v.Pos + int64(refLen) - 1
	if GenomicToCDS(delStartGenomic, t) > 0 || GenomicToCDS(delEndGenomic, t) > 0 {
		return 0, 0, 0 // touches coding sequence; the exonic path handles it
	}

	delSeq := v.Ref
	next := delEndGenomic + 1
	if t.IsReverseStrand() {
		delSeq = ReverseComplement(delSeq)
		next = delStartGenomic - 1
	}
	exonStart := GenomicToCDS(next, t)
	if exonStart < 1 {
		return 0, 0, 0
	}

	k := 0
	for {
		idx := int(exonStart) - 1 + k
		if idx >= len(t.CDSSequence) || t.CDSSequence[idx] != delSeq[k%refLen] {
			break
		}
		k++
	}
	if k < refLen {
		return 0, 0, 0
	}
	cdsStart = exonStart + int64(k-refLen)
	return cdsStart, cdsStart + int64(refLen) - 1, k
}

// cdsPosRangeStr builds a CDS position string like "c.34del", "c.34_36del",
// "c.35dup", etc. using a stack-allocated buffer for zero intermediate allocations.
func cdsPosRangeStr(start, end int, suffix string) string {
	var buf [64]byte
	n := copy(buf[:], "c.")
	n += putInt64(buf[n:], int64(start))
	if start != end {
		buf[n] = '_'
		n++
		n += putInt64(buf[n:], int64(end))
	}
	n += copy(buf[n:], suffix)
	return string(buf[:n])
}

// shiftInsertionBuf shifts an insertion rightward (3' direction) in CDS space.
// It operates on the provided mutable byte slice, avoiding allocations.
// The shift stops at maxIdx (0-based) to avoid crossing exon-exon junctions.
// Returns the new anchor index.
func shiftInsertionBuf(seq []byte, cdsAnchorIdx int, cdsSeq string, maxIdx int) int {
	idx := cdsAnchorIdx
	for idx < maxIdx && idx+1 < len(cdsSeq) && cdsSeq[idx+1] == seq[0] {
		first := seq[0]
		copy(seq, seq[1:])
		seq[len(seq)-1] = first
		idx++
	}
	return idx
}

// genomicToHGVScPos converts a genomic position to an HGVSc position string.
// For coding transcripts, returns strings like "76", "88+1", "89-2", "-14", "*6".
// For non-coding transcripts, returns transcript-relative positions like "127", "42+5".
func genomicToHGVScPos(pos int64, t *cache.Transcript) string {
	// Check if position is in an exon
	exon := t.FindExon(pos)
	if exon != nil {
		if t.IsProteinCoding() {
			return exonicHGVScPos(pos, exon, t)
		}
		return nonCodingExonicPos(pos, t)
	}

	// Intronic position — find flanking exons
	return intronicHGVScPos(pos, t)
}

// nonCodingExonicPos returns the transcript position for an exonic position in a
// non-coding transcript. Position 1 is the first exonic base at the 5' end.
func nonCodingExonicPos(pos int64, t *cache.Transcript) string {
	txPos := GenomicToTranscriptPos(pos, t)
	if txPos > 0 {
		return strconv.FormatInt(txPos, 10)
	}
	return ""
}

// exonicHGVScPos returns the HGVSc position for an exonic position.
func exonicHGVScPos(pos int64, exon *cache.Exon, t *cache.Transcript) string {
	// Check if in CDS
	if exon.IsCoding() && pos >= exon.CDSStart && pos <= exon.CDSEnd {
		cdsPos := GenomicToCDS(pos, t)
		if cdsPos > 0 {
			return strconv.FormatInt(cdsPos, 10)
		}
	}

	// 5'UTR or 3'UTR
	if t.IsForwardStrand() {
		if pos < t.CDSStart {
			return fiveprimeUTRPos(pos, t)
		}
		if pos > t.CDSEnd {
			return threeprimeUTRPos(pos, t)
		}
	} else {
		if pos > t.CDSEnd {
			return fiveprimeUTRPos(pos, t)
		}
		if pos < t.CDSStart {
			return threeprimeUTRPos(pos, t)
		}
	}

	return ""
}

// fiveprimeUTRPos returns the HGVSc position for a 5'UTR position (e.g., "-14").
// Counts exonic bases from the position to CDS start.
func fiveprimeUTRPos(pos int64, t *cache.Transcript) string {
	dist := exonicDistance(pos, t.CDSStart, t) // forward strand
	if t.IsReverseStrand() {
		dist = exonicDistance(t.CDSEnd, pos, t) // reverse strand: CDSEnd is 5' end
	}
	if dist < 0 {
		dist = -dist
	}
	return "-" + strconv.FormatInt(dist, 10)
}

// threeprimeUTRPos returns the HGVSc position for a 3'UTR position (e.g., "*6").
// Counts exonic bases from CDS end to the position.
func threeprimeUTRPos(pos int64, t *cache.Transcript) string {
	dist := exonicDistance(t.CDSEnd, pos, t) // forward strand
	if t.IsReverseStrand() {
		dist = exonicDistance(pos, t.CDSStart, t)
	}
	if dist < 0 {
		dist = -dist
	}
	return "*" + strconv.FormatInt(dist, 10)
}

// exonicDistance counts the number of exonic bases between two genomic positions.
// Both positions are inclusive. Only counts bases that fall within exons.
// from < to in genomic coordinates.
func exonicDistance(from, to int64, t *cache.Transcript) int64 {
	if from > to {
		from, to = to, from
	}
	var dist int64
	for _, exon := range t.Exons {
		// Overlap between [from, to] and [exon.Start, exon.End]
		overlapStart := from
		if exon.Start > overlapStart {
			overlapStart = exon.Start
		}
		overlapEnd := to
		if exon.End < overlapEnd {
			overlapEnd = exon.End
		}
		if overlapStart <= overlapEnd {
			dist += overlapEnd - overlapStart + 1
		}
	}
	// Subtract 1 because we want distance not span (from position is the anchor)
	if dist > 0 {
		dist--
	}
	return dist
}

// intronicHGVScPos computes the HGVSc position for an intronic position.
// Format: c.{boundary}+{offset} or c.{boundary}-{offset}
// Uses FindNearestExonIdx for O(log n) flanking exon lookup.
func intronicHGVScPos(pos int64, t *cache.Transcript) string {
	idx := t.FindNearestExonIdx(pos)
	n := len(t.Exons)
	if n == 0 {
		return ""
	}

	// Find flanking exons by checking idx and its neighbors.
	// Exons are sorted ascending by Start. An intronic position has
	// an upstream exon (End < pos) and a downstream exon (Start > pos).
	var upstreamExon, downstreamExon *cache.Exon
	for _, i := range [3]int{idx - 1, idx, idx + 1} {
		if i < 0 || i >= n {
			continue
		}
		exon := &t.Exons[i]
		if exon.End < pos {
			if upstreamExon == nil || exon.End > upstreamExon.End {
				upstreamExon = exon
			}
		}
		if exon.Start > pos {
			if downstreamExon == nil || exon.Start < downstreamExon.Start {
				downstreamExon = exon
			}
		}
	}

	if upstreamExon == nil && downstreamExon == nil {
		return ""
	}

	// Calculate distances to each flanking exon boundary
	var distToUpstream, distToDownstream int64
	if upstreamExon != nil {
		distToUpstream = pos - upstreamExon.End
	}
	if downstreamExon != nil {
		distToDownstream = downstreamExon.Start - pos
	}

	// Pick the closer exon boundary
	useUpstream := true
	if upstreamExon == nil {
		useUpstream = false
	} else if downstreamExon != nil && distToDownstream < distToUpstream {
		useUpstream = false
	}

	if t.IsForwardStrand() {
		if useUpstream {
			boundaryPos := exonBoundaryHGVScPos(upstreamExon.End, upstreamExon, t)
			return boundaryPos + "+" + strconv.FormatInt(distToUpstream, 10)
		}
		boundaryPos := exonBoundaryHGVScPos(downstreamExon.Start, downstreamExon, t)
		return boundaryPos + "-" + strconv.FormatInt(distToDownstream, 10)
	}

	// Reverse strand
	if useUpstream {
		boundaryPos := exonBoundaryHGVScPos(upstreamExon.End, upstreamExon, t)
		return boundaryPos + "-" + strconv.FormatInt(distToUpstream, 10)
	}
	boundaryPos := exonBoundaryHGVScPos(downstreamExon.Start, downstreamExon, t)
	return boundaryPos + "+" + strconv.FormatInt(distToDownstream, 10)
}

// exonBoundaryHGVScPos returns the HGVSc position string for an exon boundary.
// The boundary is the exon position closest to the intron.
func exonBoundaryHGVScPos(genomicPos int64, exon *cache.Exon, t *cache.Transcript) string {
	// Non-coding transcripts use transcript-relative position
	if !t.IsProteinCoding() {
		return nonCodingExonicPos(genomicPos, t)
	}

	// If the boundary is in the CDS, use CDS position
	if exon.IsCoding() && genomicPos >= exon.CDSStart && genomicPos <= exon.CDSEnd {
		cdsPos := GenomicToCDS(genomicPos, t)
		if cdsPos > 0 {
			return strconv.FormatInt(cdsPos, 10)
		}
	}

	// If in 5'UTR region
	if t.IsForwardStrand() {
		if genomicPos < t.CDSStart {
			return fiveprimeUTRPos(genomicPos, t)
		}
		if genomicPos > t.CDSEnd {
			return threeprimeUTRPos(genomicPos, t)
		}
	} else {
		if genomicPos > t.CDSEnd {
			return fiveprimeUTRPos(genomicPos, t)
		}
		if genomicPos < t.CDSStart {
			return threeprimeUTRPos(genomicPos, t)
		}
	}

	return ""
}

// cdsExonEndIdx returns the 0-based CDS index of the last base in the same exon
// as the given 0-based CDS index. This is used to prevent 3' shifting from crossing
// exon-exon junctions per HGVS rules. Returns len(cdsSeq)-1 if no exon boundary
// information is available.
func cdsExonEndIdx(cdsIdx0 int, t *cache.Transcript) int {
	if len(t.CDSRegions) == 0 {
		return len(t.CDSSequence) - 1
	}
	for _, r := range t.CDSRegions {
		regionLen := int(r.GenomicEnd - r.GenomicStart + 1)
		regionStart := int(r.CDSOffset)
		regionEnd := regionStart + regionLen - 1
		if cdsIdx0 >= regionStart && cdsIdx0 <= regionEnd {
			return regionEnd
		}
	}
	return len(t.CDSSequence) - 1
}

// shiftInsertionThreePrime shifts an insertion rightward (3' direction) in CDS space.
// cdsAnchorIdx is the 0-based CDS index of the VCF anchor base (insertion is after this position).
// maxIdx limits the shift to avoid crossing exon-exon junctions.
// Returns the shifted inserted sequence and new anchor index.
func shiftInsertionThreePrime(insertedSeq string, cdsAnchorIdx int, cdsSeq string, maxIdx int) (string, int) {
	seq := []byte(insertedSeq)
	idx := cdsAnchorIdx
	for idx < maxIdx && idx+1 < len(cdsSeq) && cdsSeq[idx+1] == seq[0] {
		first := seq[0]
		copy(seq, seq[1:])
		seq[len(seq)-1] = first
		idx++
	}
	return string(seq), idx
}

// shiftDeletionThreePrime shifts a deletion rightward (3' direction) in CDS space.
// delStart and delEnd are 0-based inclusive CDS indices of the deleted bases.
// The shift stops at maxIdx (0-based) to avoid crossing exon-exon junctions.
// Returns the shifted start and end indices.
func shiftDeletionThreePrime(delStart, delEnd int, cdsSeq string, maxIdx int) (int, int) {
	for delEnd < maxIdx && delEnd+1 < len(cdsSeq) && cdsSeq[delStart] == cdsSeq[delEnd+1] {
		delStart++
		delEnd++
	}
	return delStart, delEnd
}

// dupResult holds the CDS position range of a detected duplication.
type dupResult struct {
	isDup    bool
	cdsStart int64 // 1-based CDS position of first duplicated base
	cdsEnd   int64 // 1-based CDS position of last duplicated base
}

// checkSpliceJunctionDup detects duplications at exon-intron boundaries where
// the VCF anchor is intronic (CDS=0) but the adjacent exonic base matches the
// inserted sequence. This handles two patterns:
//   - Forward strand c.X-1_XinsN: anchor (v.Pos) intronic, v.Pos+1 exonic
//   - Reverse strand c.X_X+1insN: anchor (v.Pos) exonic, v.Pos+1 intronic
//
// In both cases, the CDS-based dup detection in formatHGVScInsertion skips
// these because GenomicToCDS returns 0 for the intronic position.
func checkSpliceJunctionDup(v *vcf.Variant, t *cache.Transcript, insertedSeq string) dupResult {
	if len(insertedSeq) == 0 || len(t.CDSSequence) == 0 {
		return dupResult{}
	}

	// Try both flanking positions to find which one is exonic
	cdsLeft := GenomicToCDS(v.Pos, t)
	cdsRight := GenomicToCDS(v.Pos+1, t)

	// We only handle the case where exactly one side is exonic
	if (cdsLeft > 0) == (cdsRight > 0) {
		return dupResult{}
	}

	seqLen := len(insertedSeq)

	if cdsRight > 0 {
		// Forward strand pattern: anchor is intronic, next position is exonic.
		// The exonic boundary is v.Pos+1, mapping to cdsRight.
		// Check if CDS bases starting at cdsRight match the inserted sequence.
		startIdx := int(cdsRight - 1) // 0-based
		endIdx := startIdx + seqLen
		if endIdx <= len(t.CDSSequence) && t.CDSSequence[startIdx:endIdx] == insertedSeq {
			return dupResult{
				isDup:    true,
				cdsStart: cdsRight,
				cdsEnd:   cdsRight + int64(seqLen) - 1,
			}
		}
	} else {
		// Reverse strand pattern: anchor is exonic, next position is intronic.
		// The exonic boundary is v.Pos, mapping to cdsLeft.
		// Check if CDS bases ending at cdsLeft match the inserted sequence.
		endIdx := int(cdsLeft) // exclusive, since cdsLeft is 1-based and we want [start, cdsLeft]
		startIdx := endIdx - seqLen
		if startIdx >= 0 && t.CDSSequence[startIdx:endIdx] == insertedSeq {
			return dupResult{
				isDup:    true,
				cdsStart: int64(startIdx + 1),
				cdsEnd:   cdsLeft,
			}
		}
	}

	return dupResult{}
}

// checkDuplication checks if an insertion duplicates adjacent reference bases.
// Per HGVS convention, the inserted sequence must match either the bases
// immediately before or immediately after the insertion point.
// Returns the CDS positions of the duplicated bases for direct HGVS formatting.
func checkDuplication(v *vcf.Variant, t *cache.Transcript, insertedSeq string) dupResult {
	if len(insertedSeq) == 0 || len(t.CDSSequence) == 0 {
		return dupResult{}
	}

	cdsPos := GenomicToCDS(v.Pos, t)
	if cdsPos < 1 {
		return dupResult{}
	}

	cdsIdx := int(cdsPos - 1) // 0-based index of the VCF anchor base
	seqLen := len(insertedSeq)

	// Check if inserted bases match the bases at the anchor position
	// (the anchor and preceding bases in CDS).
	startIdx := cdsIdx - seqLen + 1
	if startIdx >= 0 && cdsIdx+1 <= len(t.CDSSequence) {
		if t.CDSSequence[startIdx:cdsIdx+1] == insertedSeq {
			return dupResult{
				isDup:    true,
				cdsStart: int64(startIdx + 1), // 1-based
				cdsEnd:   cdsPos,
			}
		}
	}

	// Check if inserted bases match the bases immediately after the anchor
	// (the bases right after the insertion point in CDS).
	afterStart := cdsIdx + 1
	afterEnd := afterStart + seqLen
	if afterStart >= 0 && afterEnd <= len(t.CDSSequence) {
		if t.CDSSequence[afterStart:afterEnd] == insertedSeq {
			return dupResult{
				isDup:    true,
				cdsStart: int64(afterStart + 1), // 1-based
				cdsEnd:   int64(afterEnd),
			}
		}
	}

	return dupResult{}
}
