#!/usr/bin/env python3
"""
Bidirectional nucleotide hash function
by Cathal Garvey. Code is AGPL, concept is public domain, but kindly do
attribute and let me know if you find it useful.
Will probably be rewritten in a more efficient language sometime soon.
"""

def chunks(li, ch):
    for i in range(0, len(li), ch):
        yield li[i:i+ch]

def xorfold(num, bits):
    """
    Takes a number of bit length bits and outputs one of half that bitlength
    by xoring the upper half with the lower. This is used to avoid the bias
    towards higher-value numbers introduced by the orientation canonicalisation
    function.
    """
    assert bits % 2 == 0, "Input bits cannot be an odd number."
    return (num>>bits) ^ (num & ((1<<bits)-1))

def reverse_complement(seq):
    revseq = ''.join({"A":"T","C":"G","G":"C","T":"A"}[n] for n in seq[::-1])
    return revseq

def canonical_orientation_dna(seq, values={'A':0,'C':1,'G':2,'T':3}):
    """
    Take a DNA string, and calculate its reverse complement. Based on the value
    of the initial nucleotides in either case, choose either original or complement.
    """
    revseq = reverse_complement(seq)
    for n,p in enumerate(zip(seq,revseq)):
        i,j = p
        if n > 1+(len(seq)//2):
            break
        if values[i] > values[j]:
            return seq
        elif values[j] > values[i]:
            return revseq
        else:
            continue
    # Palindrome
    return seq


def numerify_DNA(seq, inputbits, values={'A':0,'C':1,'G':2,'T':3}):
    """
    Canonicalise either forward or reverse using canonical_orientation_dna.
    Then, create a number by valuing and bit-shifting according to position,
    with values A=0,C=1,G=2,T=3 and bit-shifts of 2 positions each. So, a
    four-character string becomes an 8-bit number.
    The result should be a number with bit length equal to twice nucleotide length,
    which is equal in value for a sequence and its reverse complement.
    So: GTCG is processed as follows:
        GTCG | CGAC : CGAC
        Values      : 1201 (12 is a better lede than 10 so chose complement)
        Bitshift by : 6420
        Results     : 64 + 32 + 0 + 1 = 97
    Finally, xor the upper half of the number with the lower to get the
    final value; this avoids biasing the number due to the canonicalisation
    protocol of selecting "higher" value directions. It also halves bit length
    in the output.
    """
    seq = canonical_orientation_dna(seq)
    num = 0
    for i,n in enumerate(reversed(seq)):
        num += values[n]<<i
    xornum = xorfold(num, inputbits//2)
    return xornum

def hashDNA2(DNAseq, offset=0, bitlen = 64):
    """
    bitlen represents the desired output bit size.
    Because bit length of encoded DNA is twice the letter length of the
    pre-encoding DNA, this means that bitlen*2 is the chunking length of
    the DNA hashing window.
    However, due to the canonicalisation function which picks a DNA orientation
    to operate upon, there is an innate bias towards higher-value numbers,
    and to account for this the upper half of each encoded chunk is xored with
    the lower to get the final value, so we return to the original bitlen.
    """
    #assert len(DNAseq) % bitlen == 0, "For testing we're assuming DNA lengths evenly factorisable by half bitlen."
    accum = (2**bitlen) - 1  # Full bitmap of bitlen. Should it be 0 instead?
    # Hash offset block, if any:
    accum ^= numerify_DNA(DNAseq[:offset], bitlen*2)
    for block in chunks(DNAseq[offset:], bitlen*2):
        # Order-sensitive.
        block_value = numerify_DNA(block, bitlen*2)
        # Order-insensitive.
        accum ^= block_value
    return accum

def hashDNA(DNAseq, seq_length, bitlen = 64):
    """
    This generates two simultaneous numbers of bitlen-bits, one starting
    at zero bytes and the other at an offset so that it "closes" on the
    final byte of DNAseq. These are xored to generate the final output.
    """
    offset = seq_length % bitlen
    hash1 = hashDNA2(DNAseq, 0, bitlen)
    hash2 = hashDNA2(DNAseq, offset, bitlen)
    return hash1 ^ hash2

# == Test code ==
if __name__ == "__main__":
    import random
    test_bitlength = 64
    
    def randomDNA(l):
        return "".join((random.choice("ACGT") for i in range(l)))

    hashDNA = hashDNA2
    seen = {}
    for i in range(1000):
        testcase = randomDNA(random.randint(150,2000))
        h = hashDNA(testcase, len(testcase), test_bitlength)
        h2 = hashDNA(reverse_complement(testcase), len(testcase), test_bitlength)
        assert h == h2, "Hashes don't match in reverse for n={}: {} vs {} - for sequence: \n{}".format(i, h, h2, testcase)
        if h in seen:
            print("Collision detected for hash: ", h, "on iteration", i)
            print("Colliding sequences:")
            print("\t", testcase)
            print("\t", seen[h])
            break
        seen[h] = testcase
