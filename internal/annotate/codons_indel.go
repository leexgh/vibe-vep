package annotate

// formatCodonChangeCDS renders the Codons string Ensembl VEP reports for an
// indel or multi-nucleotide variant, working in CDS coordinates.
//
// VEP expands the affected range out to whole-codon boundaries, then writes
// ref/alt over that window with the bases that actually changed in uppercase
// and the untouched flanking bases in lowercase. Unlike a substitution the two
// sides differ in length:
//
//	c.5946del          agT/ag           (one base deleted from the codon end)
//	c.1856del          cCt/ct           (deleted base mid-codon)
//	c.4049_4051del     gAGGcc/gcc       (three bases spanning two codons)
//	c.29_30insCAG      cag/caCAGg       (insertion inside a codon)
//	c.627delinsAT      atG/atAT         (one base replaced by two)
//	c.1132_1133delinsT AGt/Tt
//
// delStart/delEnd are 1-based inclusive CDS positions of the deleted bases; an
// empty range (delEnd < delStart) means a pure insertion after delEnd. inserted
// holds the inserted bases on the coding strand, empty for a pure deletion.
//
// An insertion that lands exactly on a codon boundary disturbs no existing
// codon, and VEP reports it as "-/INSERTED".
func formatCodonChangeCDS(cdsSeq string, delStart, delEnd int64, inserted string) string {
	n := int64(len(cdsSeq))
	if n == 0 {
		return ""
	}
	pureInsertion := delEnd < delStart
	if pureInsertion && delEnd%3 == 0 {
		// Codon-boundary insertion: no existing codon is disturbed.
		var buf [256]byte
		if len(inserted)+2 > len(buf) {
			return "-/" + upper(inserted)
		}
		k := copy(buf[:], "-/")
		k += copyUpper(buf[k:], inserted)
		return string(buf[:k])
	}

	first, last := delStart, delEnd
	if pureInsertion {
		first, last = delEnd, delEnd
	}
	if first < 1 || last > n {
		return ""
	}

	codonStart := ((first-1)/3)*3 + 1
	codonEnd := ((last-1)/3 + 1) * 3
	if codonEnd > n {
		codonEnd = n
	}
	if codonStart < 1 || codonStart > codonEnd {
		return ""
	}

	// Worst case: the whole window twice, plus the insertion twice, plus '/'.
	need := int(codonEnd-codonStart+1)*2 + len(inserted) + 1
	var stack [256]byte
	buf := stack[:]
	if need > len(buf) {
		buf = make([]byte, need)
	}
	k := 0

	// ref side: deleted bases uppercase, flanks lowercase
	for pos := codonStart; pos <= codonEnd; pos++ {
		b := cdsSeq[pos-1]
		if !pureInsertion && pos >= delStart && pos <= delEnd {
			buf[k] = upperByte(b)
		} else {
			buf[k] = lowerByte(b)
		}
		k++
	}
	buf[k] = '/'
	k++
	altStart := k

	// alt side: inserted bases uppercase, retained flanks lowercase
	for pos := codonStart; pos <= codonEnd; pos++ {
		if !pureInsertion && pos == delStart && inserted != "" {
			k += copyUpper(buf[k:], inserted)
		}
		if !pureInsertion && pos >= delStart && pos <= delEnd {
			continue
		}
		buf[k] = lowerByte(cdsSeq[pos-1])
		k++
		if pureInsertion && pos == delEnd {
			k += copyUpper(buf[k:], inserted)
		}
	}
	// A window that is wholly deleted leaves the alt side empty; VEP writes
	// that as "-", e.g. "GGC/-".
	if k == altStart {
		buf[k] = '-'
		k++
	}
	return string(buf[:k])
}

func formatCodonChangeMNV(cdsSeq string, start, end int64, altBases string) string {
	n := int64(len(cdsSeq))
	if n == 0 || start < 1 || end > n || end < start {
		return ""
	}
	if int64(len(altBases)) != end-start+1 {
		return ""
	}
	codonStart := ((start-1)/3)*3 + 1
	codonEnd := ((end-1)/3 + 1) * 3
	if codonEnd > n {
		codonEnd = n
	}
	width := int(codonEnd - codonStart + 1)
	var stack [256]byte
	buf := stack[:]
	if width*2+1 > len(buf) {
		buf = make([]byte, width*2+1)
	}
	k := 0
	for pos := codonStart; pos <= codonEnd; pos++ {
		b := cdsSeq[pos-1]
		if pos >= start && pos <= end {
			buf[k] = upperByte(b)
		} else {
			buf[k] = lowerByte(b)
		}
		k++
	}
	buf[k] = '/'
	k++
	for pos := codonStart; pos <= codonEnd; pos++ {
		if pos >= start && pos <= end {
			buf[k] = upperByte(altBases[pos-start])
		} else {
			buf[k] = lowerByte(cdsSeq[pos-1])
		}
		k++
	}
	return string(buf[:k])
}

// copyUpper writes s into dst uppercased and returns the number of bytes written.
func copyUpper(dst []byte, s string) int {
	for i := 0; i < len(s); i++ {
		dst[i] = upperByte(s[i])
	}
	return len(s)
}

func upperByte(b byte) byte { return b &^ 0x20 }
func lowerByte(b byte) byte { return b | 0x20 }

func upper(s string) string {
	b := []byte(s)
	for i := range b {
		b[i] = upperByte(b[i])
	}
	return string(b)
}
