# DNA2Way
## A hash function for DNA that generates same output for forward-bias or reverse complement.
by Cathal Garvey; implementation is AGPL, concept is public domain. 
Please let me know if this is useful, and I'd appreciate attribution.

### Contains
* A python module with simple tests when run as a script
* A Go library with tests.

### Concept
Hashing functions are handy for dealing with data; they let you see when
you've seen a piece of data before without requiring you to compare painstakingly
against everything you've got on record. By compressing data down to a small
size, you can check it against a much smaller list of snippets, and treat these
as representing the larger dataset. This is one use of a "hash function";
the generation of small "hashes" that try to represent a larger and more
meaningful piece of data in a smaller amount of space, while preserving
as much uniqueness as possible.

In bioinformatics, we deal with oodles of data, and one of the idiosyncracies
of that data is that it can be represented backwards and inverted, but still
be the "same" data. That is, DNA is a bi-directional code, it has "forward"
and "reverse complement" strands (which is called which is often subjective)
with different codes, but in the real world they are just two halves of the
same molecule.

Betimes, it would be useful to hash DNA in such a way that the hashing function
could produce the same output if it encountered the reverse complement of something
later on. And, it'd be nice to do this without requiring you to separately
hash both directions of a molecule and store them together, because sometimes
we're streaming data.

This is an imperfect attempt to solve that problem, and it's only a few hours
old, so it surely has loads of room for improvement.

### Design
#### Goals:

* Our function must generate constant-size output that's reasonably well-distributed
  across bit-space for a range of possible inputs.
* Our function must be able to stream over a sequence without holding it all in
  memory, seeking backwards, or requiring a double-take in both directions.
* Our function must generate different results for different molecules, even if
  they have the same gross base composition (on average all molecules are 0.25 A,C,G,T
  so this ought to go without saying!).
* Our function must generate the same output for a molecule whether in the forward or
  reverse direction.

#### Limitations accepted in this version:

* This function currently expects that the length of the data to be hashed is known;
  this information is necessary to offset one of the streaming "windows" for digesting
  the data correctly at the outset.
* The final function should be streaming, but to get things done quickly this
  generates two hashes sequentially and xors the outputs. So it does actually read
  sequences twice, but in the same direction merely at offsets from one another,
  so it could easily be re-written using a Python generator pattern (or a file-like
  reading pattern in another language) to generate two simultaneous hashes in
  a stream-like manner.

#### Future directions:

* In extremis, a function without knowledge of the final length of the sequence
  could generate many simultaneous digests according to the desired bit-length,
  and select the correct offset at the end when the length is finally known. This
  would increase memory and CPU consumption of the hash function itself, but would
  permit hashing of data without knowing data length in advance, which might have
  some advantages for hashing huge or non-locally-held data.
  
### Function Operation
1. Data is read into the function using a "window" that is twice the desired bit length.
2. In each window, the DNA is represented as a number. The number is encoded by selecting
   either the forward or reverse bias of that windowed fragment of DNA according to which
   orientation generates the "higher" value number. To eliminate bias, the number's upper
   bits are then xored with its lower bits, generating a number of the originally desired
   bit length.
3. Each chunk of the desired bit-length is xored with an accumulating number, which begins
   as a full bitfield for the bit-length desired. Xor is an order-insensitive operation, but
   by generating large windowed chunks and xor-folding each chunk it is expected that a
   long enough stream of unique numbers is generated in the sub-operation of the function
   to create unique output for nearly any sequence of sufficient size (modulo bit-length of course).
4. To account for varying lengths of sequence, which will rarely be of the precise
   product-of-double-desired-bit-length required to evenly generate output, two digests as steps 1-3
   are conducted; one at offset 0 in the sequence, and another at an offset that means the final
   windowed chunk closes on the final byte (this offset window consumes an incomplete chunk first,
   then complete chunks thereafter, whereas the first window consumes complete chunks until the
   last chunk).
5. The two numbers generated simultaneously in (4) are xored together to get the output.
   
### Usage
Jeese, I dunno. Ask a real bioinformatician.
