use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use rust_htslib::bam::{IndexedReader, Read, Header};
use rust_htslib::bam;
use std::collections::HashMap;
use strand_specifier_lib::{Strand, LibType, check_flag};
use std::str::FromStr;
use std::cmp;
use rust_htslib::bam::record::Record;
use CigarParser::cigar::Cigar;

/// Extracts sequence names from a BAM file header.
///
/// Reads the BAM file header and returns a vector of all sequence names (chromosome names)
/// defined in the header's SN (sequence name) fields.
///
/// # Arguments
///
/// * `bam_path` - Path to the BAM file
///
/// # Returns
///
/// A `PyResult` containing a vector of sequence names as strings
///
/// # Example
///
/// ```python
/// sequences = get_header("alignment.bam")
/// # Returns: ["chr1", "chr2", "chr3", ...]
/// ```
#[pyfunction]
fn get_header(bam_path: String) -> PyResult<Vec<String>>{

    let bam = bam::Reader::from_path(&bam_path).unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut seq: Vec<String> = Vec::new(); 

    for (key, records) in header.to_hashmap() {
            for record in records {
                if record.contains_key("SN"){
                    seq.push(record["SN"].to_string());
                }
            }
    }
    Ok(seq)
}


/// Calculates coverage at each position using the pileup algorithm.
///
/// Retrieves primary aligned reads from a BAM file within a specified genomic region
/// and calculates per-base coverage, taking into account strand specificity and
/// library type. Uses the pileup approach to count reads at each position.
///
/// # Arguments
///
/// * `start` - Start position of the region (0-based, inclusive)
/// * `end` - End position of the region (0-based, exclusive)
/// * `chrom` - Chromosome/sequence name
/// * `strand` - Strand specification ("Plus", "Minus", or "NA" for unstranded)
/// * `bam_path` - Path to the indexed BAM file
/// * `lib` - Library type specification (e.g., "RF", "FR", "unstranded")
/// * `mapq_thr` - Minimum mapping quality threshold
///
/// # Returns
///
/// A `PyResult` containing a vector of coverage values (u32) for each position in the region
///
/// # Notes
///
/// - Filters out deletions, reference skips, and secondary alignments (flag 256)
/// - Applies mapping quality filtering
/// - Handles strand-specific coverage based on library type
/// #[deprecated(since="0.5.0", note="please use `cover_from_intervall` instead")]
#[pyfunction]
fn get_coverage(start:i64, end:i64, chrom: String, strand: String,
     bam_path: String, lib: String, mapq_thr: u8) -> PyResult<Vec<u32>>{
    let mut container = vec![0; (end - start) as usize];
    let mut bam = IndexedReader::from_path(&bam_path).unwrap();

    let lib_type = LibType::from(lib.as_str());
    let strand_feature = Strand::from(strand.as_str());

    bam.fetch((&chrom, start, end)).unwrap();
    let mut read_strand: Strand = Strand::Plus; 
    let mut cpt = 0;
    for p in bam.pileup() {
        let pileup = p.unwrap();
        cpt = 0;
        if start <= i64::from(pileup.pos()) && i64::from(pileup.pos()) < end {
            for alignment in pileup.alignments() {
                if !alignment.is_del() && !alignment.is_refskip() && !check_flag(alignment.record().flags(), 256, 0) && !(alignment.record().mapq() < mapq_thr){
                    
                    if strand_feature == Strand::NA{
                        cpt += 1;
                        continue
                    }

                    if let Some(read_strand) = lib_type.get_strand(alignment.record().flags()){
                        if read_strand == Strand::NA {
                            cpt += 1;
                            continue
                        }
                        else if strand_feature == read_strand{
                            cpt += 1;
                            continue
                        }
                    }
                }
            }

            container[(pileup.pos() as i64 - start as i64) as usize ] = cpt;
            
        }

    } 

    return Ok(container)

}





/// Calculates coverage at each position using an interval-based algorithm.
///
/// Alternative coverage calculation method that processes complete read alignments
/// rather than using pileup. Computes coverage by determining which positions
/// each read covers based on its CIGAR string, with support for custom flag filtering
/// and strand-specific counting.
///
/// # Arguments
///
/// * `start` - Start position of the region (0-based, inclusive)
/// * `end` - End position of the region (0-based, exclusive)
/// * `chrom` - Chromosome/sequence name
/// * `strand` - Strand specification ("Plus", "Minus", or "NA" for unstranded)
/// * `bam_path` - Path to the indexed BAM file
/// * `lib` - Library type specification (e.g., "RF", "FR", "unstranded")
/// * `mapq_thr` - Minimum mapping quality threshold (0 to disable filtering)
/// * `flag_in` - SAM flags that must be present (bitwise AND)
/// * `flag_exclude` - SAM flags that must be absent (bitwise AND)
///
/// # Returns
///
/// A `PyResult` containing a vector of coverage values (u32) for each position in the region
///
/// # Notes
///
/// - Uses CIGAR string parsing to determine read coverage intervals
/// - More flexible flag filtering compared to `get_coverage`
/// - May be more efficient for sparse coverage regions
#[pyfunction]
fn get_coverage_algo2(start:i64, end:i64, chrom: String, strand: String,
     bam_path: String, lib: String, mapq_thr: u8, flag_in: u16, flag_exclude: u16) -> PyResult<Vec<u32>>{
    ///
    let mut container = vec![0; (end - start) as usize];
    let mut bam = IndexedReader::from_path(&bam_path).unwrap();

    let lib_type = LibType::from(lib.as_str());
    let strand_feature = Strand::from(strand.as_str());

    bam.fetch((&chrom, start, end)).unwrap();
    let mut read_strand: Strand = Strand::Plus; 
    let mut cpt = 0;

    let mut record: Record;
    let mut pos_s: i64;
    let mut pos_e: i64;
    let mut cig: Cigar;
    let mut flag: u16;

    for p in bam.records() {
        record = p.unwrap();

        pos_s = record.pos();
        cig = Cigar::from_str(&record.cigar().to_string()).unwrap();
        pos_e = cig.get_end_of_aln(&pos_s);
        flag = record.flags();
        if check_flag(flag, flag_in, flag_exclude) && (mapq_thr == 0 || !(record.mapq() < mapq_thr)){
            match lib_type{
                LibType::Unstranded | LibType::PairedUnstranded => {
                    cover_from_intervall(&mut container, start, end, cig.get_reference_cover(pos_s));
                },
                _ => {
                    if let Some(read_strand) = lib_type.get_strand(flag){
                        if strand_feature == read_strand{
                            cover_from_intervall(&mut container, start, end, cig.get_reference_cover(pos_s));
                            }
                        }
                    }
            }
            //if let Some(read_strand) = lib_type.get_strand(flag){
            //    if strand_feature == read_strand{
            //        cover_from_intervall(&mut container, start, end, cig.get_reference_cover(pos_s));
            //    }
            //}
        }
    } 
    return Ok(container)
}


/// Updates a coverage vector with read coverage intervals.
///
/// Helper function that increments coverage counts for positions that overlap
/// between a feature region and read coverage intervals. Used internally by
/// coverage calculation algorithms.
///
/// # Arguments
///
/// * `feature_cover` - Mutable reference to the coverage vector to update
/// * `cover_start` - Start position of the feature region
/// * `cover_end` - End position of the feature region
/// * `read_cover` - Vector of alternating start/end positions representing read coverage intervals
///
/// # Notes
///
/// - `read_cover` format: [start1, end1, start2, end2, ...]
/// - Only positions within both the feature region and read intervals are incremented
/// - Uses efficient clamping to handle partial overlaps
pub fn cover_from_intervall(feature_cover : &mut Vec<u32>,
    cover_start: i64,
    cover_end: i64, 
    read_cover: Vec<i64>) -> (){
    
    let mut end: i64 = 0;
    
    for (i, start) in read_cover.iter().enumerate().step_by(2){
         end = read_cover[i + 1];
         if !((*start > cover_end) | (end < cover_start)){
            for ii in (*cmp::max(start, &cover_start) as usize)..(*cmp::min(&end, &cover_end) as usize){
                 feature_cover[ii - cover_start as usize] += 1;
             } 
             
         }
    }
()
}




/// A Python module for calculating BAM coverage using Rust.
///
/// This module provides efficient functions for computing per-base coverage from BAM files
/// with support for strand-specific libraries and flexible filtering options.
///
/// # Functions
///
/// - `get_header`: Extract sequence names from BAM header
/// - `get_coverage`: Calculate coverage using pileup algorithm
/// - `get_coverage_algo2`: Calculate coverage using interval-based algorithm
#[pymodule]
fn Rust_covpyo3(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(get_header, m)?)?;
    m.add_function(wrap_pyfunction!(get_coverage_algo2, m)?)?;
    Ok(())
}
