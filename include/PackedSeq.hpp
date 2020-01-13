/*
 * 
 * Copyright (c) 2019, Ritu Kundu and Joshua Casey
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** Class PackedSequence.
 * It represents a sequence with a DNA base with 4-bits letter.
 * Either 2-bits or 4-bits packing
 * 4-bits Valid letters: A,C, G, T,N
 * 2-bits valid letter
 */
#pragma once
#ifndef PACKED_SEQ_HPP
#define PACKED_SEQ_HPP
#include "globalDefs.hpp"

#define MAX_LEN_LIMIT 0xffffffffu


namespace hypo
{
const std::string c2BStrings[256] = {
    "AAAA", "AAAC", "AAAG", "AAAT",  "AACA", "AACC", "AACG", "AACT",  "AAGA", "AAGC", "AAGG", "AAGT",  "AATA", "AATC", "AATG", "AATT", 
	"ACAA", "ACAC", "ACAG", "ACAT",  "ACCA", "ACCC", "ACCG", "ACCT",  "ACGA", "ACGC", "ACGG", "ACGT",  "ACTA", "ACTC", "ACTG", "ACTT", 
	"AGAA", "AGAC", "AGAG", "AGAT",  "AGCA", "AGCC", "AGCG", "AGCT",  "AGGA", "AGGC", "AGGG", "AGGT",  "AGTA", "AGTC", "AGTG", "AGTT", 
	"ATAA", "ATAC", "ATAG", "ATAT",  "ATCA", "ATCC", "ATCG", "ATCT",  "ATGA", "ATGC", "ATGG", "ATGT",  "ATTA", "ATTC", "ATTG", "ATTT", 
	"CAAA", "CAAC", "CAAG", "CAAT",  "CACA", "CACC", "CACG", "CACT",  "CAGA", "CAGC", "CAGG", "CAGT",  "CATA", "CATC", "CATG", "CATT", 
	"CCAA", "CCAC", "CCAG", "CCAT",  "CCCA", "CCCC", "CCCG", "CCCT",  "CCGA", "CCGC", "CCGG", "CCGT",  "CCTA", "CCTC", "CCTG", "CCTT", 
	"CGAA", "CGAC", "CGAG", "CGAT",  "CGCA", "CGCC", "CGCG", "CGCT",  "CGGA", "CGGC", "CGGG", "CGGT",  "CGTA", "CGTC", "CGTG", "CGTT", 
	"CTAA", "CTAC", "CTAG", "CTAT",  "CTCA", "CTCC", "CTCG", "CTCT",  "CTGA", "CTGC", "CTGG", "CTGT",  "CTTA", "CTTC", "CTTG", "CTTT", 
	"GAAA", "GAAC", "GAAG", "GAAT",  "GACA", "GACC", "GACG", "GACT",  "GAGA", "GAGC", "GAGG", "GAGT",  "GATA", "GATC", "GATG", "GATT", 
	"GCAA", "GCAC", "GCAG", "GCAT",  "GCCA", "GCCC", "GCCG", "GCCT",  "GCGA", "GCGC", "GCGG", "GCGT",  "GCTA", "GCTC", "GCTG", "GCTT", 
	"GGAA", "GGAC", "GGAG", "GGAT",  "GGCA", "GGCC", "GGCG", "GGCT",  "GGGA", "GGGC", "GGGG", "GGGT",  "GGTA", "GGTC", "GGTG", "GGTT", 
	"GTAA", "GTAC", "GTAG", "GTAT",  "GTCA", "GTCC", "GTCG", "GTCT",  "GTGA", "GTGC", "GTGG", "GTGT",  "GTTA", "GTTC", "GTTG", "GTTT", 
	"TAAA", "TAAC", "TAAG", "TAAT",  "TACA", "TACC", "TACG", "TACT",  "TAGA", "TAGC", "TAGG", "TAGT",  "TATA", "TATC", "TATG", "TATT", 
	"TCAA", "TCAC", "TCAG", "TCAT",  "TCCA", "TCCC", "TCCG", "TCCT",  "TCGA", "TCGC", "TCGG", "TCGT",  "TCTA", "TCTC", "TCTG", "TCTT", 
	"TGAA", "TGAC", "TGAG", "TGAT",  "TGCA", "TGCC", "TGCG", "TGCT",  "TGGA", "TGGC", "TGGG", "TGGT",  "TGTA", "TGTC", "TGTG", "TGTT", 
	"TTAA", "TTAC", "TTAG", "TTAT",  "TTCA", "TTCC", "TTCG", "TTCT",  "TTGA", "TTGC", "TTGG", "TTGT",  "TTTA", "TTTC", "TTTG", "TTTT"

};

const std::string c4BStrings[256] = {
    "AA", "AC", "AG", "AT",  "AN", "AN", "AN", "AN",  "AN", "AN", "AN", "AN",  "AN", "AN", "AN", "AN", 
	"CA", "CC", "CG", "CT",  "CN", "CN", "CN", "CN",  "CN", "CN", "CN", "CN",  "CN", "CN", "CN", "CN", 
	"GA", "GC", "GG", "GT",  "GN", "GN", "GN", "GN",  "GN", "GN", "GN", "GN",  "GN", "GN", "GN", "GN", 
	"TA", "TC", "TG", "TT",  "TN", "TN", "TN", "TN",  "TN", "TN", "TN", "TN",  "TN", "TN", "TN", "TN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN", 
	"NA", "NC", "NG", "NT",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN",  "NN", "NN", "NN", "NN"
};

// Struct representing Homopolymer Compressed String
 using HPCString = struct HPCString{
	 std::string hpc;
	 std::vector<UINT32> counts;
 };

template <int NB>
class PackedSeq
{
public:
    static const BYTE byte_mask; 
	PackedSeq(): _data(),_remainder(0), _valid(true){}
    PackedSeq(const std::string &data);

	// Creates packed seq from kseq.h reading format.
    PackedSeq(const kseq_t * ks);

	// Creates packed seq from htslib seq format.
	PackedSeq(const UINT32 seq_len, const UINT32 offset, const UINT8 * hts_seq); // makes packed seq from htslib seq format.

	// Creates new packed seq from left_ind(inclusive) to right_ind(exclusive) of another packed seq (of the same type/NB)
	PackedSeq(const PackedSeq & ps, const size_t left_ind, const size_t right_ind);

	// Creates new packed seq (2bits base) from left_ind(inclusive) to right_ind(exclusive) of another packed seq (4-bits base)
	template<int MB>
	PackedSeq(const PackedSeq<MB> & ps, const size_t left_ind, const size_t right_ind) ;
	
	PackedSeq(const PackedSeq &) = default;
    PackedSeq &operator=(const PackedSeq &) = delete;
	PackedSeq(PackedSeq&&) = default;
    PackedSeq& operator=(PackedSeq&&) = default;
    ~PackedSeq() = default;

 	// Needed to convert PackedSeq<4> to PackedSeq<2>
	friend class PackedSeq<2>;

	inline bool is_valid() {return _valid;}
	// Unpacks the whole sequence
    inline std::string unpack() const {return unpack(0,get_seq_size());}
	// Unpacks the sequence from left_ind(inclusive) to right_ind(exclusive)
	std::string unpack(const size_t left_ind, const size_t right_ind) const;
	// Mutiply (data_size-1) bytes to number of bases in a byte; Add remainder
	// * 4 (=> << 2) and * 2 (=> << 1)
	inline size_t get_seq_size() const {return _data.size()>0?(((_data.size()-1)  << _ind_shift) + _remainder) : 0;}
	// Get relevant byte by dividing ind by 4 (if NB=2) or 2 (if NB=4) 
	// Get the relevant base in byte by ind % 4 (if NB=2) or ind%2 (if NB=4)
	// Shift the relevant base-bits to the end of byte.
	// And with byte_mask to keep only the relevant bits
	inline BYTE enc_base_at(size_t ind) const {assert(ind < get_seq_size()); return BYTE(_data[ind >> _ind_shift] >> (_bit_shift[ind & _mod_mask])) &  byte_mask;}
	inline BYTE base_at(size_t ind) const {return cCode[enc_base_at(ind)];}

	/** Finds the first/last position of the given kmer in the given range, if any
	 *  Range: left_ind(inclusive) to right_ind(exclusive)
	 * 	First if is_first is set; otherwise last
	 *  Stores index in result, if found
	 *  Returns whether found or not
	 * */
	bool find_kmer(const UINT64 target_kmer, const UINT k, const size_t left_ind, const size_t right_ind, const bool is_first, size_t& result) const;

	/** Checks whether the kmer at the given position is same as the given kmer
	 *  Returns whether same or not
	 * */
	bool check_kmer(const UINT64 target_kmer, const UINT k, const size_t ind)const;

	/** Finds the first/last position of the given minimser kmer (canonical) in the given range, if any
	 *  Range: left_ind(inclusive) to right_ind(exclusive)
	 * 	First if is_first is set; otherwise last
	 *  Stores index in result, if found
	 *  Returns whether found or not
	 * */
	bool find_canonical_kmer(const UINT64 target_canonical_kmer, const UINT k, const size_t left_ind, const size_t right_ind, const bool is_first, size_t& result) const;

	/** Checks whether the minimiser kmer at the given position is same as the given minimiser kmer (canonical)
	 *  Returns whether same or not
	 * */
	bool check_canonical_kmer(const UINT64 target_canonical_kmer, const UINT k, const size_t ind)const;

private:
    std::vector<BYTE> _data;
    UINT8 _remainder; // Number of bases in the last byte
	bool _valid;
	static const UINT8 _bases_in_byte;
	static const UINT8 _ind_shift; // Find byte number by /2(equivalent to >> 1) if NB is 4 or /4(equivalent to >>2) if NB is 2
	static const size_t _mod_mask; // %4(for NB=2) and %2(For NB=4) is equi to extracting the last 2 or 1 bits/bit resp.
	static const UINT8 _bit_shift[]; // Shift the bits by this amount to get the relevent base at the end.
	static const UINT8 _cHts_char[16]; // Shift the bits by this amount to get the relevent base at the end.
}; // class PackedSeq

// TODO: Make it elegant. Currently, this is just a work-around to avoid error "specialization of ... after instantiation ..."

template <> template<>
PackedSeq<2>::PackedSeq(const PackedSeq<4> & ps, const size_t left_ind, const size_t right_ind);

extern template class PackedSeq<2>;
extern template class PackedSeq<4>;

} // namespace hypo
#endif