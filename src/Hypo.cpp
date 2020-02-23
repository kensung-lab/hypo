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
/** Defines the class Hypo.
 * It is the master class carrying out the polishing.
 */


#include "Hypo.hpp"
#include "suk/SolidKmers.hpp"
#include "Window.hpp"
# include <omp.h>



namespace hypo
{
    
Hypo::Hypo(const InputFlags& flags): _cFlags(flags), _monitor() {
    omp_set_num_threads(_cFlags.threads);
}
void Hypo::polish() {
    std::ofstream gStagefile(STAGEFILE, std::ofstream::out | std::ofstream::app);
    if (_cFlags.intermed && !gStagefile.is_open()) {
        fprintf(stderr, "[Hypo::Hypo] Error: File open error: Stage File (%s) exists but could not be opened!\n",STAGEFILE);
        exit(1);
    }
    auto num_threads = _cFlags.threads;
    
    ///////////////////////////////////
    /* Get Solid kmers */
    std::unique_ptr<suk::SolidKmers> pSK;
    _monitor.start();
    if (_cFlags.done_stage< STAGE_SK) {
        pSK = std::make_unique<suk::SolidKmers>(_cFlags.k);
        bool is_success = pSK->initialise(_cFlags.sr_filenames,num_threads, _cFlags.sz_in_gb, _cFlags.cov, true, AUX_DIR);
        if (!is_success) {
            fprintf(stderr, "[Hypo::SolidKmers] Error: KMC Output: Could not have successful run of SUK for computing Solid kmers!\n");
            exit(1);
        }
        if (_cFlags.intermed) {
            if (pSK->store(SKFILE)) {
                std::string tm = _monitor.stop("[Hypo:Hypo]: Computed Solid kmers. ");
                gStagefile << "Stage:SolidKmers [" << tm << "]\t" << STAGE_SK << std::endl;
            }
            else {
                fprintf(stderr, "[Hypo::SolidKmers] Error: File Saving: Could not store the DS for Solid kmers!\n");
                exit(1);
            }
        }
        else {
            _monitor.stop("[Hypo:Hypo]: Computed Solid kmers. ");
        }
    }
    else {
        pSK = std::make_unique<suk::SolidKmers>(_cFlags.k);
        if (!(pSK->load(SKFILE))) {
                fprintf(stderr, "[Hypo::SolidKmers] Error: File Loading: Could not load the DS for Solid kmers (%s)!\n",SKFILE);
                exit(1);
            }
        _monitor.stop("[Hypo:Hypo]: Loaded Solid kmers. ");
    }
    fprintf(stdout, "[Hypo::Hypo] Info: Number of (canonical) solid kmers (nonhp) : %lu\n",pSK->get_num_solid_kmers()); 
    
    ///////////////////////////////////
    /* Get contigs */
    _monitor.start();
    UINT32 cid = 0;
    gzFile fp = gzopen(_cFlags.draft_filename.c_str(), "r");
    kseq_t *seq;
    seq = kseq_init(fp);
    int l;
    while((l = kseq_read(seq)) >= 0) {
        _cname_to_id[std::string(seq->name.s,seq->name.l)] = cid;
        _contigs.emplace_back(std::make_unique<Contig>(cid, seq));
        ++cid;
    }
    _monitor.stop("[Hypo:Hypo]: Loaded Contigs. ");
    _contigs.shrink_to_fit();
    _alignment_store.resize(_contigs.size());
    ///////////////////////////////////
    /* Find solid positions in contigs */
    // TODO: Try parallelising it at segment level
    _monitor.start();
    #pragma omp parallel for schedule(static,1) 
    for (UINT64 i = 0; i < _contigs.size(); ++i) {        
        _contigs[i]->find_solid_pos(pSK);
    }
     _monitor.stop("[Hypo:Hypo]: Found Solid pos in contigs. ");

    /** /////////////////////// BATCH WISE PROCESSING BEGINS //////////////////////////// **/
    _contig_batch_size = (_cFlags.processing_batch_size==0)?(_contigs.size()):(_cFlags.processing_batch_size);
    UINT32 num_batches = _contigs.size()/_contig_batch_size;
    if (_contigs.size()%_contig_batch_size!=0) {
        ++num_batches;
    }
    fprintf(stdout, "[Hypo::Hypo] Info: Number.of contigs: %lu; Number of batches: %u\n",_contigs.size(),num_batches);
    _sf_short = sam_open(_cFlags.sr_bam_filename.c_str(), "r");
    _sam_header_short = sam_hdr_read(_sf_short);
    _hts_align_short = bam_init1();

    if (_cFlags.lr_bam_filename!="") {
        _sf_long = sam_open(_cFlags.lr_bam_filename.c_str(), "r");
        _sam_header_long = sam_hdr_read(_sf_long);
        _hts_align_long = bam_init1();
    }
    for (UINT32 batch_id = 0; batch_id< num_batches; ++batch_id) {
        fprintf(stdout, "********** [Hypo::Hypo] Info: BATCH-ID: %u\n",batch_id);
        UINT32 initial_cid = batch_id * _contig_batch_size;
        UINT32 final_cid = std::min(UINT32(_contigs.size()),initial_cid + _contig_batch_size);
        ////////////////////////////)///////
        /* Get alignments */
        _monitor.start();
        create_alignments(true, batch_id);
        _monitor.stop("[Hypo:Hypo]: Loaded alignments. ");

        ///////////////////////////////////
        /* Find strong regions (SR) in contigs and divide WR into windows */
        _monitor.start();
        for (UINT64 cid=initial_cid; cid < final_cid; ++cid) {
            if (cid==initial_cid || cid==final_cid-1 || cid%10==0) {
                std::cout << "Kmer support update: " << cid <<std::endl;
            }
                
            #pragma omp parallel for
            for (UINT64 t = 0; t < _alignment_store[cid].size(); ++t) { 
                _alignment_store[cid][t]->update_solidkmers_support(_cFlags.k, *(_contigs[cid]));
            }            
        }  
        _monitor.stop("[Hypo:Hypo]: Solid kmers support update. ");

        _monitor.start();
        #pragma omp parallel for schedule(static,1) 
        for (UINT64 i = initial_cid; i < final_cid; ++i) {
            _contigs[i]->prepare_for_division(_cFlags.k);
        }
        UINT64 num_sr =0;
        UINT64 len_sr = 0;
        for (UINT64 i = initial_cid; i < final_cid; ++i) {
            num_sr += _contigs[i]->get_num_sr();
            len_sr += _contigs[i]->get_len_sr();
        }
        fprintf(stdout, "[Hypo::Hypo] Info: Total number of SR: %lu; Total length of SR: %lu\n",num_sr,len_sr); 
        _monitor.stop("[Hypo:Hypo]: Finding SR (and preparing for division). ");  

        _monitor.start();
        for (UINT64 cid=initial_cid; cid < final_cid; ++cid) {
            if (cid==initial_cid || cid==final_cid-1 || cid%10==0) {
                std::cout << "Minimser support update: " << cid <<std::endl;
            }
            #pragma omp parallel for
            for (UINT64 t = 0; t < _alignment_store[cid].size(); ++t) { 
                _alignment_store[cid][t]->update_minimisers_support(*(_contigs[cid]));
            }            
        }  
        _monitor.stop("[Hypo:Hypo]: Minimisers support update. ");
        
        _monitor.start();
        #pragma omp parallel for schedule(static,1) 
        for (UINT64 i = initial_cid; i < final_cid; ++i) {
            _contigs[i]->divide_into_regions();
        }
        _monitor.stop("[Hypo:Hypo]: Division into windows. ");

        ///////////////////////////////////
        /* Fill windows with short arms */
        _monitor.start();
        for (UINT64 cid=initial_cid; cid < final_cid; ++cid) {
            if (cid==initial_cid || cid==final_cid-1 || cid%10==0) {
                std::cout << "SA: " << cid <<std::endl;
            }
            #pragma omp parallel for
            for (UINT64 i = 0; i < _alignment_store[cid].size(); ++i) { 
                _alignment_store[cid][i]->find_short_arms(_cFlags.k, *(_contigs[cid]));
            }
        }
        _monitor.stop("[Hypo:Hypo]: Short arms computing. ");
        _monitor.start();
        // This will destroy alignments after use.
        #pragma omp parallel for schedule(static,1)  
        for (UINT64 i = initial_cid; i < final_cid; ++i) {
                    _contigs[i]->fill_short_windows(_alignment_store[i]);
                    _alignment_store[i].clear();
        }
        _monitor.stop("[Hypo:Hypo]: Short arms filling. ");
        ///////////////////////////////////
        /* Get long alignments */
        if (_cFlags.lr_bam_filename!="") {
            _monitor.start();
            create_alignments(false,batch_id);
            _monitor.stop("[Hypo:Hypo]: Loaded alignments of Long reads. ");

            ///////////////////////////////////
            /* Get long arms */
            _monitor.start();
            #pragma omp parallel for schedule(static,1) 
            for (UINT64 i = initial_cid; i < final_cid; ++i) {
                _contigs[i]->prepare_long_windows();
            }
            for (UINT64 cid=initial_cid; cid < final_cid; ++cid) {
                #pragma omp parallel for
                for (UINT64 i = 0; i < _alignment_store[cid].size(); ++i) { 
                    _alignment_store[cid][i]->find_long_arms(*(_contigs[cid]));
                }
            }

            // This will destroy alignments after use.
            #pragma omp parallel for schedule(static,1)  
            for (UINT64 i = initial_cid; i < final_cid; ++i) {
                _contigs[i]->fill_long_windows(_alignment_store[i]);
                _alignment_store[i].clear();
            }
            _monitor.stop("[Hypo:Hypo]: Long arms filling. ");
        }
        else {
            Contig::set_no_long_reads();
        }         
        
        ///////////////////////////////////
        /* Polish with arms */
        _monitor.start();
        Window::prepare_for_poa(_cFlags.score_params,_cFlags.threads);
        for (UINT64 i = initial_cid; i < final_cid; ++i) {
            UINT64 num_reg = _contigs[i]->get_num_regions();
            #pragma omp parallel for schedule(static,1) 
            for (UINT64 w=0; w < num_reg; ++w) {
                if (_contigs[i]->is_valid_window(w)) {
                    auto tid = omp_get_thread_num();
                    _contigs[i]->generate_consensus(w,tid);     
                }
            }
        }
        _monitor.stop("[Hypo:Hypo]: POA of windows. ");
    }
    /** /////////////////////// BATCH WISE PROCESSING BEGINS //////////////////////////// **/
     _alignment_store.clear();
    _alignment_store.shrink_to_fit(); 

    ///////////////////////////////////
    /* Write results */
    _monitor.start();
    std::ofstream ofile(_cFlags.output_filename);
    if (!ofile.is_open()) {
        fprintf(stderr, "[Hypo::Hypo] Error: File open error: Output File (%s) could not be opened!\n",_cFlags.output_filename.c_str());
        exit(1);
    }
    //std::ofstream bedfile(BEDFILE);
    for (UINT64 i = 0; i < _contigs.size(); ++i) {
       ofile << *(_contigs[i]);
       //_contigs[i]->generate_inspect_file(bedfile);
    }
    ofile.close();
    _monitor.stop("[Hypo:Hypo]: Writing results. ");

    gStagefile.close();
    //bedfile.close();
    _monitor.total("[Hypo:Hypo]: Overall. ");

    // Clean-up
    _contigs.clear();
}

void Hypo::create_alignments(bool is_sr, UINT32 batch_id) {
    UINT8 mq = _cFlags.map_qual_th;

    auto sf = (is_sr) ? (_sf_short) : (_sf_long);
    auto sam_header = (is_sr) ? (_sam_header_short) : (_sam_header_long);
    auto hts_align = (is_sr) ? (_hts_align_short) : (_hts_align_long);
    
    UINT32 initial_cid = batch_id * _contig_batch_size;
    UINT32 final_cid = initial_cid + _contig_batch_size;

    UINT64 num_invalid = 0;
    UINT64 num_alns = 0;
    UINT64 num_aln_read = 0;
    UINT64 progress= 0;
    while(sam_read1(sf, sam_header, hts_align)>=0) {
        if ((num_aln_read%100000000)==0) {
            std::cout << "Aln processed (in 100 M): " << progress <<std::endl;
            ++progress;
        }
        ++num_aln_read;
        // Ignore if unmapped/secondary/duplicate mapping or with failed QC
        if (hts_align->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {continue;}
        // Ignore if mapping quality <=1
        if (hts_align->core.qual < mq) {continue;}
        std::string cname(sam_hdr_tid2name(sam_header,hts_align->core.tid));
        if (_cname_to_id.find(cname) == _cname_to_id.end()) {
            fprintf(stderr, "[Hypo::Hypo] Error: Alignment File error: Contig-reference (%s) does not exist in the draft!\n",cname.c_str());
            exit(1);
        }
        UINT32 cid = _cname_to_id[cname];
        if (is_sr) {  // short reads     
            _alignment_store[cid].emplace_back(std::make_unique<Alignment>(*(_contigs[cid]), hts_align));
        }
        else { // long reads
            _alignment_store[cid].emplace_back(std::make_unique<Alignment>(*(_contigs[cid]), _cFlags.norm_edit_th, hts_align));
        }
        if (!((_alignment_store[cid].back())->is_valid)) {
            _alignment_store[cid].pop_back();
            --num_alns;
            ++num_invalid;
        }
        ++num_alns;
        if (cid >= final_cid) { // cid belonging to next batch begins
            break;
        }        
    }
    fprintf(stdout, "[Hypo::Hypo] Info: Number of alignments (Batch %u): loaded (%lu) invalid (%lu)\n",batch_id,num_alns,num_invalid); 
    if (final_cid>=_contigs.size()) { // last batch done; clean up
        bam_destroy1(hts_align);
        sam_close(sf);
    }   
}

} // namespace hypo
