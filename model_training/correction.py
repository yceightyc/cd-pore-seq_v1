import os
import re
import itertools
import pysam
import parasail
import sys
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ProcessPoolExecutor, as_completed

input_fasta = sys.argv[1]
input_bam = sys.argv[2]
output_folder = sys.argv[3]
batch_size = 10000

os.makedirs(output_folder, exist_ok=True)

# Get bam
existing_batches = {
    int(fname.split("_")[1].split(".")[0])  
    for fname in os.listdir(output_folder)
    if fname.startswith("batch_") and fname.endswith(".bam")
}

# Convert CIGAR string to pysam format
def cigarstring_to_pysam(cigar_str):
    op_map = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    return [(op_map[m.group(2)], int(m.group(1)))
            for m in re.finditer(r'(\d+)([MIDNSHP=X])', cigar_str)]

# Get Fasta records
def batch_iterator(iterator, batch_size):
    while True:
        batch = list(itertools.islice(iterator, batch_size))
        if not batch:
            break
        yield batch

# Load bam into dictionary
def load_bam_reads(bam_path):
    read_lookup = {}
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for read in bamfile.fetch(until_eof=True):
            read_lookup[read.query_name] = read.query_sequence
    return read_lookup

# Alignment correction
def process_batch(args):
    batch_index, batch_records, read_lookup, output_folder = args
    read_data = []
    ref_lengths = {}

    for record in batch_records:
        full_id = record.id
        ref_seq = str(record.seq)

        if "_" not in full_id:
            continue

        contig, query_name = full_id.split("_", 1)
        ref_name = full_id

        read_seq = read_lookup.get(query_name)
        if not read_seq:
            continue

        if read_seq == ref_seq:
            cigar_tuples = [(0, len(read_seq))]
        else:
            alignment = parasail.nw_trace_scan(read_seq, ref_seq, 10, 1, parasail.dnafull)
            cigar_str = alignment.cigar.decode.decode("utf-8")
            cigar_tuples = cigarstring_to_pysam(cigar_str)

        read_data.append({
            'query_name': query_name,
            'query_seq': read_seq,
            'ref_name': ref_name,
            'ref_seq': ref_seq,
            'cigar': cigar_tuples
        })

        if ref_name not in ref_lengths:
            ref_lengths[ref_name] = len(ref_seq)

    if not read_data:
        print(f"Batch {batch_index:03d} had no matches — skipped.")
        return

    # Build header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': name, 'LN': ln} for name, ln in ref_lengths.items()]
    }

    output_bam = os.path.join(output_folder, f"batch_{batch_index:03d}.bam")
    with pysam.AlignmentFile(output_bam, "wb", header=header) as bamfile:
        for record in read_data:
            a = pysam.AlignedSegment()
            a.query_name = record['query_name']
            a.query_sequence = record['query_seq']
            a.flag = 0
            a.reference_id = bamfile.get_tid(record['ref_name'])
            a.reference_start = 0
            a.mapping_quality = 60
            a.cigar = record['cigar']
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            a.query_qualities = pysam.qualitystring_to_array("I" * len(record['query_seq']))
            bamfile.write(a)

    print(f"Wrote {len(read_data)} alignments to {output_bam}")

# Parallelisation
if __name__ == "__main__":
    print("Preloading reads from BAM...")
    read_lookup = load_bam_reads(input_bam)
    print(f"Loaded {len(read_lookup)} reads")

    print("Beginning batch processing...")
    print(f"Found {len(existing_batches)} completed batches. Resuming...")

    fasta_iter = SeqIO.parse(input_fasta, "fasta")
    all_batches = []

    for i, batch in enumerate(batch_iterator(fasta_iter, batch_size), start=1):
        if i in existing_batches:
            continue  # skip already done
        all_batches.append((i, batch, read_lookup, output_folder))

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_batch, args) for args in all_batches]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print("ERROR:", e)
