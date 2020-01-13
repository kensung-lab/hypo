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
/** Defines the class PackedSequence.
 * It represents a sequence with a DNA base with 4-bits letter.
 * Either 2-bits or 4-bits packing
 * 4-bits Valid letters: A,C, G, T,N
 * 2-bits valid letter
 */
#include "PackedSeq.hpp"

namespace hypo
{

template <int NB>
const UINT8 PackedSeq<NB>::byte_mask = BYTE((NB==2)?(0x03):(0x0f));

template <int NB>
const UINT8 PackedSeq<NB>::_bases_in_byte = UINT8((NB==2)?(4):(2));

// Find byte number by /2(equivalent to >> 1) if NB is 4 or /4(equivalent to >>2) if NB is 2
// Find numbases by *2(equivalent to << 1) if NB is 4 or *4(equivalent to <<2) if NB is 2
template <int NB>
const UINT8 PackedSeq<NB>::_ind_shift = UINT8((NB==2)?(2):(1));

// Shift the bits by this amount to get the relevent base at the end.
//template <int NB>
//const std::vector<UINT8> PackedSeq<NB>::_bit_shift = (NB==2)?{6,4,2,0}:{4,0}; 
template <>
const UINT8 PackedSeq<2>::_bit_shift[] = {6,4,2,0};

template <>
const UINT8 PackedSeq<4>::_bit_shift[] = {4,0};

// %4(for NB=2) and %2(For NB=4) is equivalent to extracting the last 2 or 1 bits/bit resp
// & with this mask to get the mod
template <int NB>
const size_t PackedSeq<NB>::_mod_mask = size_t((NB==2)?(3):(1));

template <int NB>
const UINT8 PackedSeq<NB>::_cHts_char[16] = {4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4};

template <int NB>
PackedSeq<NB>::PackedSeq(const std::string& data):_valid(true) {
    static_assert(NB==2 || NB==4,"[Hypo::PackedSequence] Packed sequence base can only be 2 or 4.");
    if (data.size() > MAX_LEN_LIMIT) {
        fprintf(stderr, "[Hypo::PackedSeq] Error: Length exceed limit: The length of a sequence is %lu which exceeds the limit of %u !\n",data.size(),MAX_LEN_LIMIT);
        exit(1);
    }

    auto seq_len = data.size();
    UINT32 num_bytes = (seq_len >> _ind_shift); //equi to divide by 4 (NBis 2) or 2 (if NB is 4)
    if ((seq_len & _mod_mask) == 0) { // completely fits
        _remainder = _bases_in_byte;
    }
    else{
        _remainder = UINT8(seq_len & _mod_mask);
        ++num_bytes;
    }    
    _data.resize(num_bytes,0);
    UINT8 shift = UINT8(8-NB);
    for (UINT32 i=0; i < seq_len; ++i) {
        BYTE b = cNt4Table[BYTE(data[i])];
        if (NB==2 && b > 3) {
            fprintf(stderr, "[Hypo::PackedSeq] Error: Wrong base (Can not pack in 2 bits): Base %c in a sequence is not A, C, G, or T !\n",data[i]);
            _valid=false;
            break;
        }
        UINT32 curr_byte = UINT32(i >> _ind_shift);
        _data[curr_byte] |= (b << shift);
        shift = (shift==0) ? (UINT8(8-NB)) : (shift-NB);
    }
    _data.shrink_to_fit();
}

template <int NB>
PackedSeq<NB>::PackedSeq(const kseq_t * ks):_valid(true) {
    static_assert(NB==2 || NB==4,"[Hypo::PackedSequence] Packed sequence base can only be 2 or 4.");
    if (ks->seq.l > MAX_LEN_LIMIT) {
        fprintf(stderr, "[Hypo::PackedSeq] Error: Length exceed limit: The length of a sequence is %lu which exceeds the limit of %u !\n",ks->seq.l,MAX_LEN_LIMIT);
        exit(1);
    }
    auto seq_len = ks->seq.l;
    UINT32 num_bytes = (seq_len>> _ind_shift); //equi to divide by 4 (NBis 2) or 2 (if NB is 4)
    if ((seq_len & _mod_mask) == 0) { // completely fits
        _remainder = _bases_in_byte;
    }
    else{
        _remainder = UINT8(seq_len & _mod_mask);
        ++num_bytes;
    }   
    _data.resize(num_bytes,0);
    UINT8 shift = UINT8(8-NB);
    for (UINT i=0; i < seq_len; ++i) {
        BYTE b = cNt4Table[BYTE(ks->seq.s[i])];
        if (NB==2 && b > 3) {
            fprintf(stderr, "[Hypo::PackedSeq] Error: Wrong base (Can not pack in 2 bits): Base %c in a sequence is not A, C, G, or T !\n",ks->seq.s[i]);
            _valid=false;
            break;
        }
         UINT32 curr_byte = UINT32(i >> _ind_shift);
        _data[curr_byte] |= (b << shift);
        shift = (shift==0) ? (UINT8(8-NB)) : (shift-NB);
    }
    _data.shrink_to_fit();
}

template <int NB>
PackedSeq<NB>::PackedSeq(const UINT32 seq_len, const UINT32 offset, const UINT8 * hts_seq): _valid(true) {
    static_assert(NB==2 || NB==4,"[Hypo::PackedSequence] Packed sequence base can only be 2 or 4.");
    if (seq_len > MAX_LEN_LIMIT) {
        fprintf(stderr, "[Hypo::PackedSeq] Error: Length exceed limit: The length of a sequence is %u which exceeds the limit of %u !\n",seq_len,MAX_LEN_LIMIT);
        exit(1);
    }
    UINT32 num_bytes = (seq_len>> _ind_shift); //equi to divide by 4 (NBis 2) or 2 (if NB is 4)
    if ((seq_len & _mod_mask) == 0) { // completely fits
        _remainder = _bases_in_byte;
    }
    else{
        _remainder = UINT8(seq_len & _mod_mask);
        ++num_bytes;
    }   
    _data.resize(num_bytes,0);
    UINT8 shift = UINT8(8-NB);
    for (UINT32 i=0; i < seq_len; ++i) {
        BYTE b = _cHts_char[BYTE(bam_seqi(hts_seq,i+offset))];
        if (NB==2 && b > 3) {
            //fprintf(stderr, "[Hypo::PackedSeq] Error: Wrong base (Can not pack in 2 bits): Base code %u in a sequence is not A, C, G, or T !\n",bam_seqi(hts_seq,i+offset));
            _valid=false;
            break;
        }
         UINT32 curr_byte = UINT32(i >> _ind_shift);
        _data[curr_byte] |= (b << shift);
        shift = (shift==0) ? (UINT8(8-NB)) : (shift-NB);
    }
    _data.shrink_to_fit();
}

// Creates new packed seq (of the same type) from left_ind(inclusive) to right_ind(exclusive) of another packed seq
template <int NB>
PackedSeq<NB>::PackedSeq(const PackedSeq<NB> & ps, const size_t left_ind, const size_t right_ind) {
    //std::cout << "Left : Right : old seq len : " << left_ind << " " << right_ind <<" " <<ps.get_seq_size()<<std::endl;
    size_t old_seq_len = ps.get_seq_size();
    assert(right_ind <= old_seq_len && left_ind<=right_ind);
    size_t seq_len = right_ind-left_ind;
    UINT32 num_bytes = (seq_len>> _ind_shift); //equi to divide by 4 (NBis 2) or 2 (if NB is 4)
    if ((seq_len & _mod_mask) == 0) { // completely fits
        _remainder = _bases_in_byte;
    }
    else{
        _remainder = UINT8(seq_len & _mod_mask);
        ++num_bytes;
    }
    _data.reserve(num_bytes);
    auto start_byte_ind = left_ind >> _ind_shift;
    auto end_byte_ind = (right_ind-1) >> _ind_shift;
    auto start_ind_in_byte = left_ind & _mod_mask;
    auto last_ind_in_byte = (right_ind-1) & _mod_mask;

    if (start_ind_in_byte == 0) { // simply copy
        _data.assign(ps._data.begin()+start_byte_ind,ps._data.begin()+end_byte_ind+1);
    }
    else {
        UINT8 mult = UINT8((NB==2) ? (1) : (2));
        UINT8 lshift = UINT8 (start_ind_in_byte <<mult);  // i*2 (NB is 2) or i*4 (NB is 4)
        UINT8 rshift = UINT8 (8 - lshift); 
        BYTE rmask = BYTE((0x1 << lshift)-1);
        size_t ind = 0;
        for (size_t i=start_byte_ind; i < end_byte_ind; ++i,++ind) {
            UINT8 b = UINT8((ps._data[i] << lshift) | ((ps._data[i+1] >> rshift) & rmask));
            _data.emplace_back(b); 
        }
        // handle last byte; if last_ind_in_byte<start_ind_in_byte then already covered in th loop
        if (last_ind_in_byte>=start_ind_in_byte) {   
            _data.emplace_back(ps._data[end_byte_ind] << lshift);    
        }
    }
    _data.shrink_to_fit();  
}

// Creates new packed seq (2bit base) from left_ind(inclusive) to right_ind(exclusive) of another packed seq (4-bits base)
template <>template<>
PackedSeq<2>::PackedSeq(const PackedSeq<4> & ps, const size_t left_ind, const size_t right_ind) {
    auto old_seq_len = ps.get_seq_size();
    assert(right_ind <= old_seq_len && left_ind<=right_ind);
    size_t seq_len = right_ind-left_ind;
    UINT32 num_bytes = (seq_len>> _ind_shift); //equi to divide by 4 (NBis 2) or 2 (if NB is 4)
    if ((seq_len & _mod_mask) == 0) { // completely fits
        _remainder = _bases_in_byte;
    }
    else{
        _remainder = UINT8(seq_len & _mod_mask);
        ++num_bytes;
    }
    _data.resize(num_bytes,0);
    UINT8 shift = UINT8(6);
    size_t ps_ind = left_ind;
    
    for (UINT32 i=0; i < seq_len; ++i,++ps_ind) {
        //std:: cout << "i ps_ind: " << i << " " << std::endl;
        BYTE b = (ps_ind & 0x01)? ps._data[ps_ind >> 1] : BYTE(ps._data[ps_ind >> 1] >> 4);
        b = b &  0x0f;
        if (b > 3) {
            fprintf(stderr, "[Hypo::PackedSeq] Error: Wrong base (Can not pack in 2 bits): Base code at %lu in a sequence is not A, C, G, or T !\n",ps_ind);
            exit(1);
        }
         UINT32 curr_byte = UINT32(i >> 2); // divide by 4
        _data[curr_byte] |= (b << shift);
        shift = (shift==0) ? (UINT8(6)) : (shift-2);
    }
    _data.shrink_to_fit();
}


// Unpacks the sequence from left_ind(inclusive) to right_ind(exclusive)
template <int NB>
std::string PackedSeq<NB>::unpack(const size_t left_ind, const size_t right_ind) const {
    auto seq_len = get_seq_size();
    assert(left_ind >=0 && right_ind <= seq_len && left_ind<=right_ind);
    if (left_ind==right_ind) {
        return std::string("");
    }
    const std::string* byte_to_chars = (NB==2) ? (c2BStrings) : (c4BStrings);
    std::string unpacked;
    size_t subseq_len = right_ind-left_ind;
    unpacked.reserve(subseq_len);
    auto start_byte_ind = left_ind >> _ind_shift;
    auto end_byte_ind = (right_ind-1) >> _ind_shift;

    if (start_byte_ind<end_byte_ind) { // first and the last bytes are not the same
        // handle first byte
        auto fb = byte_to_chars[_data[start_byte_ind]];
        auto start_ind_in_byte = left_ind & _mod_mask;
        unpacked = std::string(fb,start_ind_in_byte);
        
        for (size_t i=start_byte_ind+1; i < end_byte_ind; ++i) {
            //std::cout << std::hex << (int)_data[i] << std::endl;
            unpacked += byte_to_chars[_data[i]];
        }
    }    
    // handle last byte
    auto lb = byte_to_chars[_data[end_byte_ind]];
    auto last_ind_in_byte = (right_ind-1) & _mod_mask;
    auto start_ind_in_byte = (start_byte_ind<end_byte_ind)?(0):(left_ind & _mod_mask);
    unpacked += std::string(lb,start_ind_in_byte,last_ind_in_byte+1);
    return unpacked;
}

template <int NB>
bool PackedSeq<NB>::check_kmer(const UINT64 target_kmer, const UINT k, const size_t ind)const {
    // Assumes kmer is based on 2-bit packing
    auto right_ind = ind+k;
    assert(right_ind <= get_seq_size());
    UINT64 kmer=0;
    UINT kmer_len=0;
    const UINT64 kmask = (1ULL<<2*k) - 1;
    bool found = false;
    for (size_t i=ind; i < right_ind; ++i) {
        BYTE b = enc_base_at(i);
        if (b < 4) { //ACGT
            kmer = ((kmer << 2) | b) & kmask;
            if (kmer_len<k) { ++kmer_len;}                
        }
        else { // N; reset
            kmer_len = 0;
            kmer=0;
        }
        if (kmer_len==k && kmer==target_kmer) { //found
            found = true;
            break;
        }
    }
    return found;
}

template <int NB>
bool PackedSeq<NB>::find_kmer(const UINT64 target_kmer, const UINT k, const size_t left_ind, const size_t right_ind, const bool is_first, size_t& result)const {
    // Assumes kmer is based on 2-bit packing
    //std::cout << "L : R: S " << left_ind << " " << right_ind << " " << get_seq_size() <<std::endl;
    assert(right_ind <= get_seq_size() && left_ind<=right_ind);
    if (left_ind==right_ind) {return false;}
    UINT64 kmer=0;
    UINT kmer_len=0;
    const UINT64 kmask = (1ULL<<2*k) - 1;
    bool found = false;
    for (size_t i=left_ind; i < right_ind; ++i) {
        BYTE b = enc_base_at(i);
        if (b < 4) { //ACGT
            kmer = ((kmer << 2) | b) & kmask;
            if (kmer_len<k) { ++kmer_len;}                
        }
        else { // N; reset
            kmer_len = 0;
            kmer=0;
        }
        if (kmer_len==k && kmer==target_kmer) { //found
            result = i-k+1;
            found = true;
            if (is_first) {
                break;
            }
        }
    }
    return found;
}

template <int NB>
bool PackedSeq<NB>::check_canonical_kmer(const UINT64 target_canonical_kmer, const UINT k, const size_t ind)const {
    // Assumes kmer is based on 2-bit packing
    auto right_ind = ind+k;
    assert(right_ind <= get_seq_size());
    UINT32 shift = 2 * (k - 1);
    UINT32 mask = (1ULL<<2*k) - 1;
    UINT32 kmer[2] = {0,0};
    UINT kmer_len=0;
    bool found = false;
    for (size_t i=ind; i < right_ind; ++i) {
        BYTE b = enc_base_at(i);
        if (b < 4) { //ACGT
            ++kmer_len;
            kmer[0] = (kmer[0] << 2ull | b) & mask;           // forward k-mer
            kmer[1] = (kmer[1] >> 2ull) | (3ULL^b) << shift; // reverse k-mer
            int z = kmer[0] < kmer[1] ? 0 : 1;
            if (kmer_len>=k) { 
                if (kmer[z]==target_canonical_kmer) {
                    found = true;
                    break;
                }                
            }                
        }
        else { // N; reset
            kmer_len = 0;
        }
    }
    return found;
}

template <int NB>
bool PackedSeq<NB>::find_canonical_kmer(const UINT64 target_canonical_kmer, const UINT k, const size_t left_ind, const size_t right_ind, const bool is_first, size_t& result)const {
    // Assumes kmer is based on 2-bit packing
    //std::cout << "L : R: S " << left_ind << " " << right_ind << " " << get_seq_size() <<std::endl;
    assert(right_ind <= get_seq_size() && left_ind<=right_ind);
    if (left_ind==right_ind) {return false;}
    UINT32 shift = 2 * (k - 1);
    UINT32 mask = (1ULL<<2*k) - 1;
    UINT32 kmer[2] = {0,0};
    UINT kmer_len=0;
    bool found = false;
    for (size_t i=left_ind; i < right_ind; ++i) {
        BYTE b = enc_base_at(i);
        if (b < 4) { //ACGT
            ++kmer_len;
            kmer[0] = (kmer[0] << 2ull | b) & mask;           // forward k-mer
            kmer[1] = (kmer[1] >> 2ull) | (3ULL^b) << shift; // reverse k-mer
            int z = kmer[0] < kmer[1] ? 0 : 1;
            if (kmer_len>=k) { 
                if (kmer[z]==target_canonical_kmer) {
                    result = i-k+1;
                    found = true;
                    if (is_first) {
                        break;
                    }
                } 
            }                
        }
        else { // N; reset
            kmer_len = 0;
        }
    }
    return found;
}


template class PackedSeq<2>;
template class PackedSeq<4>;

} // namespace hypo
