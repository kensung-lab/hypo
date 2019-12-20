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
/** Class Hypo.
 * It is the master class carrying out the polishing.
 */
#pragma once
#include <unordered_map>
#include <thread>
#include <memory>

#include "slog/Monitor.hpp"

#include "globalDefs.hpp"
#include "Contig.hpp"
#include "Alignment.hpp"

namespace hypo
{
class SolidKmers;

class Hypo {
public:
    Hypo(const InputFlags&);
    ~Hypo() = default;
    void polish();
private:
    const InputFlags& _cFlags;
    UINT32 _contig_batch_size;
    samFile *_sf_short;
    bam_hdr_t *_sam_header_short;
    bam1_t *_hts_align_short;
    samFile *_sf_long;
    bam_hdr_t *_sam_header_long;
    bam1_t *_hts_align_long;
    slog::Monitor _monitor;
    
    std::unordered_map<std::string, UINT32> _cname_to_id;
    std::vector<std::unique_ptr<Contig>> _contigs;
    std::vector<std::vector<std::unique_ptr<Alignment>>> _alignment_store;

    void create_alignments(bool is_sr, UINT32 batch_id); 
}; // Hypo
} // namespace hypo