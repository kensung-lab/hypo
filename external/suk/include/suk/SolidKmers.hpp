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
/** Defines the class SolidKmers.
 * It represents the solid (unique) kmers as a bit-vector.
 * DNA alphabet has four letters/bases:"ACGT". Therefore, there are `4^k` number of distinct kmers possible. 
 * DNA alphabet can be represented by 2 bits per base: `00`, `01`, `10`, `11`. Thus each kmer can be represented in `2k` bits. 
 * If we encode each base of a kmer into bits, it represents a unique number. Consequently, we have one-to-one correspondence between kmers and numbers of `2k` bits.
 * Thus, we can show the collection of kmers appearing in data as a bit-vector of length `2k`. 
 * Each index of the vector represents a kmer; if a kmer is prsent, we set the bit at the corresponding index; otherwise, reset.
 */
#pragma once

#include <sdsl/bit_vectors.hpp>

namespace suk
{

// Types
  using UINT = unsigned int;
  using INT = int;
  using INT8 = int8_t;
  using UINT8 = uint8_t;
  using BYTE = uint8_t;
  using INT16 = int16_t;
  using UINT16 = uint16_t;
  using INT32 = int32_t;
  using UINT32 = uint32_t;
  using INT64 = int64_t;
  using UINT64 = uint64_t;

const std::vector<char> cCode = {'A','C','G','T','N'};
const std::vector<char> cRCCode = {'T','G','C','A','N'};
// A: 0, C: 1, G:2, T: 3; everything else is 4.
const BYTE cNt4Table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

using CutOffs = struct SCutOffs{
    UINT err;
    UINT mean;
    UINT upper;
    UINT lower;
  };

class SolidKmers
{
public:    
    SolidKmers(const UINT k);
    SolidKmers(const SolidKmers&) = delete;
    SolidKmers& operator=(const SolidKmers&) = delete;

    SolidKmers(SolidKmers&&) = default;
    SolidKmers& operator=(SolidKmers&&) = default;

    ~SolidKmers() = default;

    /** Initialises bit-vector of solid unique kmers.
     * filenames is the vector containing the names of the read-files
     * threads is the number of threads
     * max_memory is the maximum memory usage of KMC
     * coverage is the expected coverage of the input
     * is_fastq is whether the input is in fastq format or not
     * exclude_hp if set, ignores the solid unique kmers that have homopolymers at the starting or ending
     * tmp_directory is the temporary directory for KMC to use, and also to put KMC output in
     * Returns true if the bit-vector is set successfully.
     * */
    bool initialise(const std::vector<std::string> & filenames, const UINT32 threads, const UINT32 max_memory, const UINT32 coverage, const bool exclude_hp, const std::string tmp_directory);

    /** Loads the bitvector from a previously stored file (Assumes k is the same) .
     * Returns true if the bit-vector is loaded successfully.
     * */
    bool load(std::string infile);

    /** Stores the bitvector into a file with the given name. 
     * Returns true if the bit-vector is stored successfully.
     * */
    inline bool store(std::string outfile) const {return sdsl::store_to_file(_bv, outfile);}

    /** Stores the kmers set in the bitvector (in text format) into a file with the given name. 
     * Returns true if the kmers are written successfully.
     * */
    bool dump_txt(const std::string ofile);

    /** Checks whether a kmer (in encoded format) is solid unique kmer. 
     * kmer_id is the kmer in encoded form (2bits per base encoding).
     * Returns true if the kmer_id is set.
     * */
    inline bool is_solid(UINT64 kmer_id) const {assert(kmer_id<_bv.size());return _bv[kmer_id];}

    /** Checks whether a kmer (text/not-encoded format) is solid unique kmer. 
     * 
     * kmer is the kmer in text format. (Assumes the given kmer is of length k.)
     * Returns true if the corresponding kmer-index is set in the bit-vector.
     * */
    bool is_solid(std::string kmer) const;   

    /** Returns the number of solid (canonical) unique kmers. */ 
    inline UINT64 get_num_solid_kmers() const {return _num_Solid_kmers;}

    /* Returns the value of k with which this class has been initialised */
    inline UINT get_k() const {return _k;}

private:
    const UINT _k;
    sdsl::bit_vector _bv;
    UINT64 _num_Solid_kmers;

    CutOffs find_cutoffs(const std::vector<size_t> & histArray);

}; // class PackedSeq

} // namespace hypo
