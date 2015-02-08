package dna2way

import (
	"fmt"
	"math"
)

// For when you don't want to import fmt just to make errors!
//type serror string
//func (s *serror) Error() string {
//	return *s
//}

const (
)

var (
	nucleotides = [...]string{"A","C","G","T"}
	complementMap = map[rune]rune{
			'A' : 'T',
			'C' : 'G',
			'G' : 'C',
			'T' : 'A',
	}
	nucleotideValues = map[rune]uint{
		'A': 0,
		'C': 1,
		'G': 2,
		'T': 3,
	}
)

// Xors the lowest portion of num that is _bits_ in bitlength with the
// immediately higher portion of num that is _bits_ in bitlength, after
// bit-shifting the higher portion downwards by _bits_
func xorfold(num, bits uint) (uint64, error) {
	if bits % 2 != 0 {
		return 0, fmt.Errorf("Bits must be an even number in xorfold.")
	}
	return uint64(((num>>bits) & ((1<<bits)-1)) ^ (num & ((1<<bits)-1))), nil
}

// Does what it says on the tin
func reverseComplement(nucs string) string {
	outbits := make([]rune, len(nucs))
	for ind, char := range nucs {
		revind := len(nucs) - 1 - ind 
	    outbits[revind] = complementMap[char]
	}
	return string(outbits) //strings.Join(outbits, "")
}

// Take a DNA string, and calculate its reverse complement. Based on the value
//  of the initial nucleotides in either case, choose either original or complement.
func canonicalOrientationDNA(seq string) string {
	if len(seq) == 0 {
		return ""
	}
	revseq := reverseComplement(seq)
	max := 1+len(seq)/2
	for n := 0; n < max; n++ {
		i := rune(seq[n])
		j := rune(revseq[n])
		if nucleotideValues[i] > nucleotideValues[j] {
			return seq
		} else if nucleotideValues[i] < nucleotideValues[j] {
			return revseq
		} else {
			continue
		}
	}
	// Palindrome
	return seq
}

func numerifyDNA(seq string, outputbits uint) (uint64, error) {
	seq = canonicalOrientationDNA(seq)
	l := len(seq)
	num := uint(0)
	for i, n := range seq {
		i = l - 1 - i
		num += nucleotideValues[n] << uint(i)
	}
	return xorfold(num, outputbits)
}

func stringChunks(s string, chunklen int) []string {
	output := make([]string, 0, 1 + len(s)/chunklen)
	for i := 0; i < len(s); i += chunklen {
		sliceends := int(math.Min(float64(i+chunklen), float64(len(s))))
		output = append(output, s[i : sliceends])
	}
	return output
}

func compressDNA(DNAseq string, offset int, outputbits uint) (uint64, error) {
	accum := uint64(math.Exp2(float64(outputbits))) - 1
	block_value, err := numerifyDNA(DNAseq[:offset], outputbits)
	if err != nil {
		return 0, err
	}
	accum ^= block_value
	for _, block := range stringChunks(DNAseq[offset:], int(outputbits * 2)) {
		block_value, err := numerifyDNA(block, outputbits)
		if err != nil {
			return 0, err
		}
		accum ^= block_value
	}
	return accum, nil
}

func HashDNA(DNAseq string, outputbits uint) (uint64, error) {
	hash1, err := compressDNA(DNAseq, 0, outputbits)
	if err != nil {
		return 0, err
	}
	offset := len(DNAseq) % int(outputbits)
	if offset == 0 {
		return hash1, nil
	}
	hash2, err := compressDNA(DNAseq, offset, outputbits)
	if err != nil {
		return 0, err
	}
	return hash1 ^ hash2, nil
}
