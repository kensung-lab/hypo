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
/** Class Alignment.
 * It represents an alignment and contains functionality for kmer -support update; minimizer-support update, and arms finding
 */
#pragma once
#include <sdsl/bit_vectors.hpp>
#include "globalDefs.hpp"
#include "PackedSeq.hpp"

//#include "Contig.hpp"

namespace hypo
{
class Contig;
enum class ArmType: UINT8 {
    INTERNAL, 
    PREFIX, 
    SUFFIX,
    EMPTY
};

using Arm = struct SArm{
    const UINT32 windex;
    const PackedSeq<2> arm;
    const ArmType armtype;
    SArm (const UINT32 ind, const PackedSeq<2>& ps, UINT32 left, UINT32 right, const ArmType at): windex(ind), arm(std::move(PackedSeq<2>(ps,left,right))), armtype(at){}
    SArm (const UINT32 ind): windex(ind), arm(std::move(PackedSeq<2>())), armtype(ArmType::EMPTY){}
  };

class Alignment {
public:
    // Alignment for short read (does not consider normalised edit distance)
    Alignment(Contig& contig, bam1_t *hts_align);

    // Alignment for long read (considers normalised edit distance)
    Alignment(Contig& contig, UINT64 norm_edit_th, bam1_t *hts_align);

    Alignment(const Alignment &) = delete;
    Alignment &operator=(const Alignment &) = delete;
	Alignment(Alignment&&) = delete;
    Alignment& operator=(Alignment&&) = delete;
     ~Alignment()= default;

    // Used to pop_back invalid alignment (a long read alignment may be invalid) from the alignment_store by hypo
    bool is_valid; // Can be false only for a long read

    void update_solidkmers_support (const UINT k, Contig& pc);
    void update_minimisers_support (Contig& contig);

    // For adding to short windows
    void find_short_arms(const UINT k, Contig& contig);

    // For adding to long windows
    void find_long_arms(Contig& contig);

    void add_arms(const Contig& contig);

private:
    UINT32 _rb; // ref beginning/start pos
    UINT32 _re; // ref ending pos (points to 1 past the index of last aligned base)
    std::string _qname; 
    UINT32 _qab; // aligned query beg pos (0 based) in the sequence stored
    UINT32 _qae; // aligned query endpos (0 based) in the sequence stored (points to 1 past the index of last aligned base)
    PackedSeq<2> _apseq; // aligned seq of query in Packed_Seq<2> format
    std::vector<UINT32> _cigar; // cigar string in htslib cigar format
    UINT32 _num_cigar_op; // number of cigar operations
  
    std::vector<Arm> _arms;

    void initialise_pos(const bam1_t *hts_align);
    void copy_data(const bam1_t *hts_align);

    std::vector<UINT32> find_bp(const sdsl::bit_vector::select_1_type& Sreg_pos, const std::vector<RegionType>& reg_type, const UINT32 beg_ind, const UINT32 end_ind);
    void prepare_short_arm(const UINT k, const UINT32 windex, const UINT32 qb, const UINT32 qe, const ArmType armtype, Contig& contig);
    
}; // Alignment
} // namespace hypo