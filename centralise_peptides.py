import csv, sys, os
from collections import defaultdict
import bio_util as util

SPECIES_CACHE_FILE = 'species_cache.txt'

def centralise_peptides(in_csv_path, out_path_root, align_path=None, half_width=7, pad_chr='-', id_col=0, pep_seq_col=4, sep=','):
                        
  pep_sites = defaultdict(list)
  
  alt_ids = util.get_uniprot_alt_ids()
  sub_ids = {}
  iso_ids = {}
  
  out_path_root  = os.path.splitext(out_path_root)[0]
  out_csv_path   = out_path_root + '.csv'
  out_fasta_path = out_path_root + '.fasta'
  
  id_letter = chr(ord("A")+id_col)
  pep_seq_letter = chr(ord("A")+pep_seq_col)
  
  util.info(f'Sending output to {out_path_root}.fasta')
  util.info(f'Cetralising parameters:')
  util.info(f' .. half width; {half_width}')
  util.info(f' .. pad character; "{pad_chr}"')
  util.info(f' .. protein ID column; {id_col+1} ({id_letter})')
  util.info(f' .. peptide sequence column; {pep_seq_col+1} ({pep_seq_letter})')
  
  site_counts = defaultdict(int)
  
  with open(in_csv_path, newline='') as file_obj:
    util.info(f'Reading {in_csv_path}')
    reader = csv.reader(file_obj, delimiter=sep)
    head = next(reader)
    for row in reader:
      pep_id = row[id_col]
      pid, site, null1, null2, n_peps = pep_id.split('_')
      
      if pid in alt_ids:
        aid = alt_ids[pid] # Get primary ID
        if aid not in sub_ids:
          util.info(f' .. repacing obsolete {pid} with {aid}')
        sub_ids[aid] = pid # Primary to orig, secondary
        pid = aid
      
      site_counts[(pid, site)] += 1
      
      if site_counts[(pid, site)] > 1:
        continue
        
      seqs = row[pep_seq_col].split(';')
      seqs = [x.strip('_') for x in seqs]
      
      pep_sites[pid].append((site, seqs))
      
  pids = sorted(pep_sites)
  util.info(f'Fetching {len(pids):,} sequences')
  seq_dict = util.get_uniprot_columns(pids, ['sequence'])
  alt_sites = defaultdict(list)
  valid_sites = {}
  
  for pid in pep_sites:
    if pid not in seq_dict:
      util.warn(f' .. Sequence missing for {pid}!')
      continue
  
    seq = seq_dict[pid]
    
    for site, pep_seqs in pep_sites[pid]:
      for pep_seq in pep_seqs:
        if pep_seq in seq:
          valid_sites[(pid, site)] = pep_seq
          break
      
      else:
        alt_sites[(pid, site)] = pep_seqs
        util.info(f' .. Peptide sequence missing for {pid} {site} {pep_seq}')

  if alt_sites:
    util.info(f'Fallback isoform check')
    for pid, site in alt_sites:
      pep_seqs = alt_sites[(pid, site)]
      aids = [pid+'-1', pid+'-2']
      alt_seq_dict = util.get_uniprot_columns(aids, ['sequence'])
      
      for aid in aids:
        if aid in alt_seq_dict:
          seq = alt_seq_dict[aid]
          
          for pep_seq in pep_seqs:
            if pep_seq in seq:
              valid_sites[(pid, site)] = pep_seq
              iso_ids[pid] = aid
              util.info(f' .. Isoform sequence used for {aid}')
              break
          else:
            continue
          break    
  
  # A get alignment data
  
  align_dict = defaultdict(dict)
  align_refs = {}
  
  if align_path:
    concat_align_dict = util.read_fasta(align_path)
    species_dict = {}
    
    if os.path.exists(SPECIES_CACHE_FILE):
      with open(SPECIES_CACHE_FILE) as file_obj:
        util.info(f'Reading species cache from {SPECIES_CACHE_FILE}')
        for line in file_obj:
           pid, species = line.rstrip().split('\t')
           species_dict[pid] = species
    # Split/allocate according to original queries
    n0 = len(species_dict)
    
    pids = []
    qmap = {}
    
    for head, aseq in concat_align_dict.items():
      query, pid = head.split('|')
      qid, rid = query.split(':')
      align_refs[qid] = rid
      
      if query in align_dict:
        prev = next(iter( align_dict[query]))
        check_len = len(align_dict[query][prev])
        
        if check_len != len(aseq):
          util.critical(f'Alignment error. Alignment length changed for {query} {check_len} vs {len(aseq)}')
      
      align_dict[query][pid] = aseq
        
      if pid not in species_dict:
        pids.append(pid)
      
      if len(pids) == 100:
        sd = util.get_uniprot_columns(pids, ['organism_name'])
        species_dict.update(sd)
        util.info(f' .. fetched species info for {len(species_dict):,} sequences', line_return=True)
        pids = []
       
    if pids:
      sd = util.get_uniprot_columns(pids, ['organism_name'])
      species_dict.update(sd)
      util.info(f' .. fetched species info for {len(species_dict):,} sequences', line_return=True)
      pids = []
    
    if len(species_dict) > n0:
      with open(SPECIES_CACHE_FILE, 'w') as out_file_obj:
        for pid, species in species_dict.items():
          line = f'{pid}\t{species}\n'
          out_file_obj.write(line)
 
        util.info(f'Wrote species cache to {SPECIES_CACHE_FILE}')
    
    # Filter for unconserved gaps

    for query in align_dict:
      non_gap_counts = defaultdict(int)
      n = len(align_dict[query])
      
      for pid in align_dict[query]:
        seq = align_dict[query][pid]
        for i, aa in enumerate(seq):
          if aa != '-':
            non_gap_counts[i] += 1
      
      alen = len(seq)
      keep_cols = sorted([i for i, c in non_gap_counts.items() if c >= (n/2.0)])

      if len(keep_cols) < alen:
        for pid in align_dict[query]:
          seq = align_dict[query][pid]
          align_dict[query][pid] = ''.join([seq[i] for i in keep_cols])
     
  # Write out subsequences
  
  pad = pad_chr * half_width
  subseq_dict = {}
  n_pep = 0
  
  with open(out_csv_path, 'w') as csv_ofo:
    if align_path:
      csv_ofo.write('primary_uniprot_id,orig_uniprot_id,species,num_sites,site,pep_seq\n')
    else:
      csv_ofo.write('primary_uniprot_id,orig_uniprot_id,num_sites,site,pep_seq\n')
    
    for pid, site in valid_sites:
      if pid in seq_dict:
        used_pid = iso_ids.get(pid, pid)
        orig_pid = sub_ids.get(pid, pid)
        query = f'{orig_pid}:{pid}' # secondary : primary
        
        seq = pad + seq_dict[pid] + pad
        ns = site_counts[(pid, site)]
        
        site_pos = int(site[1:])-1 # Numbering started at 1
        i = site_pos + half_width
        pep_seq = seq[i-half_width:i+half_width+1]
        
        if not pep_seq:
          util.info(f' .. primary accession sequence too short for {pid} {site}')
          continue
        
        a = seq.index(pep_seq)
        b = a + len(pep_seq)
        if not a <= site_pos < b:
          util.info(f' .. site {pid} {site} not within MS peptide sequence')
          continue
        
        if align_path:
          adict = align_dict[query]
          
          if pid not in adict: # Should now always be true
            util.info(f' .. self alignment not present for {pid} ({query})')
            continue
          
          pseq = adict[pid]
          ref_site_pos = -1
          
          k = 0
          for j, aa in enumerate(pseq):
            if aa != '-':
              if site_pos == k:
                ref_site_pos = j
              k += 1
          
          if ref_site_pos < 0:
            util.info(f' .. site {pid} {site} not in reference sequence length {k}/{len(seq)}, site_pos {site_pos}')
            continue
          
          c = max(0, ref_site_pos-half_width)
          d = min(len(pseq), ref_site_pos+half_width+1)
          
          aids = sorted(adict)
          aids.remove(pid)
          aids.insert(0, pid)
          done = set()
           
          for aid in aids:
            aseq = adict[aid]
            asite_aa = aseq[ref_site_pos]
            
            if asite_aa not in 'STY':
              #print(f' .. homolog of {pid} : {aid} aligment site {asite_aa}{ref_site_pos} not Ser, Thr or Tyr')
              continue
              
            apep_seq = aseq[c:d] # Includes gaps
            
            if apep_seq in done:
              continue
              
            done.add(apep_seq)
            
            apep_seq = apep_seq.replace('-', pad_chr)
            
            species = species_dict[aid]
            
            asite_pos = len(aseq[:ref_site_pos].replace('-','')) # Position in unaligned homologue sequence

            csv_line = f'{aid},{aid},{species},{ns},{asite_aa}{asite_pos},{apep_seq}\n'
            csv_ofo.write(csv_line)
            subseq_dict[f'{aid}_{asite_aa}{asite_pos}___{ns}|{species.split()[0]}'] = apep_seq
            n_pep += 1
       
        else:
          csv_line = f'{used_pid},{sub_ids.get(aid, aid)},{ns},{site},{pep_seq}\n'
          csv_ofo.write(csv_line)
          subseq_dict[f'{used_pid}_{site}___{ns}'] = pep_seq
          n_pep += 1
 
  
  print(f'Written {n_pep:,} rows to {out_csv_path}')
  
  util.write_fasta(out_fasta_path, subseq_dict)
        

def main(argv=None):

  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]
  
  DEFAULT_WIDTH = 7
  DEFAULT_ID_COL = '1'
  DEFAULT_PEP_SEQ_COL = '5'
  DEFAULT_PAD_CHR = 'X'
  
  epilog = 'For further help email tstevens@mrc-lmb.cam.ac.uk'

  arg_parse = ArgumentParser(prog='centralise_peptides.py', description='Extract protein subsequences around MS identified sites of interest, e.g. phosphorylation.',
                             epilog=epilog, prefix_chars='-', add_help=True)
  
  arg_parse.add_argument(metavar='CSV_FILE', dest='i',
                         help='Input CSV file containing mass spec data, must contain columns for protein ID with site location and for peptide sub-sequence.')
  
  arg_parse.add_argument('-o', '--out-file-root', default=None, metavar='OUT_FILE_ROOT', dest="o",
                         help=f'Output file path prefix (i.e. no file extension) for writing output CSV and FASTA format files. ' \
                              'If unspecified a default will be based on the input CSV file.')
  
  arg_parse.add_argument('-a', '--align-file-path', default=None, metavar='ALIGN_PATH', dest="a",
                         help=f'Optional alignment file path, in concatenated FASTA format, as output by "align_orthologues.py".')
  
  arg_parse.add_argument('-w', metavar='HALF_WIDTH', default=DEFAULT_WIDTH, type=int,
                         help=f'Number of amino acid residues to add to both sides the central site to make each subsequence. Default: {DEFAULT_WIDTH}')
  
  arg_parse.add_argument('-p', metavar='PAD_CHAR', default=DEFAULT_PAD_CHR,
                         help=f'Padding/extension character to pad the start/end of protein sequences. Default: {DEFAULT_PAD_CHR}')
  
  arg_parse.add_argument('-ci', '--id-column', default=DEFAULT_ID_COL, metavar='COLUMN', dest="ci",
                         help='Column in the CSV file containing UniProt protein IDs and site of interest in the format {PROT_ID}_{AA}{RES_NUM}___{COUNT}, e.g. "Q9Z2D1_S631___1". ' +\
                              f'Column identifier may be a number, starting from 1, of a letter starting from A. Default {DEFAULT_ID_COL}')
  
  arg_parse.add_argument('-cp', '--pep-seq-column', default=DEFAULT_PEP_SEQ_COL, metavar='COLUMN',dest="cp",
                         help=f'Column in the CSV file containing MS identified peptide subsequence(s). ' \
                              f'Column identifier may be a number, starting from 1, of a letter starting from A. Default {DEFAULT_PEP_SEQ_COL}')

  args = vars(arg_parse.parse_args(argv))

  in_csv_path   = args['i']
  out_path_root = args['o']
  align_path    = args['a'] 
  pad_char      = args['p']
  half_width    = max(0, args['w'])
  prot_id_col   = args['ci']
  pep_seq_col   = args['cp']
  
  if not out_path_root:
    if align_path:
      out_path_root = os.path.splitext(in_csv_path)[0] + f'_centralised_ortho_subseqs_w{half_width}' 
    
    else:
      out_path_root = os.path.splitext(in_csv_path)[0] + f'_centralised_subseqs_w{half_width}' 
    
  if not os.path.exists(in_csv_path):
    print(f'ERROR: File {in_csv_path} not found')
    sys.exit(1)     
  
  if prot_id_col[0].isdigit():
    prot_id_col = int(prot_id_col)-1
  else:
    prot_id_col = ord(prot_id_col[0].upper()) - ord('A')

  if pep_seq_col[0].isdigit():
    pep_seq_col = int(pep_seq_col)-1
  else:
    pep_seq_col = ord(pep_seq_col[0].upper()) - ord('A')

  if prot_id_col == pep_seq_col:
    print(f'ERROR: Protein ID and peptide subsequence columns cannot be the same')
    sys.exit(1)     
  
  centralise_peptides(in_csv_path, out_path_root, align_path, half_width, pad_char, prot_id_col, pep_seq_col)
  

if __name__ == "__main__":
  main()



