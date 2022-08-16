from collections import defaultdict
import bio_util as util

def get_prot_sequences(csv_path, id_col=0, subseq_col=3):
  
  alt_ids = util.get_uniprot_alt_ids()
  
  #def head_to_id(head):
  #  return head.split()[0].split('|')[1]
  
  ref_map = {}
  orig_ids = {}
  n_alt = 0
  
  phospho_sites = defaultdict(list)
  with open(csv_path) as file_obj:
    head = file_obj.readline()
    
    for line in file_obj:
      site_id, prot_names, gene_name, sub_seq, *null = line.split(',')
      pid, phospho_res, null1, null2, site_count = site_id.split('_')
      
      if pid in alt_ids:
        n_alt += 1
        rid = alt_ids[pid]
        orig_ids[rid] = pid
        pid = rid
        
      phospho_sites[pid].append((phospho_res, sub_seq))
  
  qids = sorted(phospho_sites)    
  pids = set()
  seq_dict = {}
  
  for pid, seq in util.get_uniprot_columns(qids, ['sequence']).items():
    key = f'{orig_ids.get(pid, pid)}:{pid}'
    seq_dict[key] = seq
    pids.add(pid)
   
  missing_ids = [pid for pid in phospho_sites if pid not in pids]
  
  # Use secondary ID mapping file, as used in CA DNN db construction
  
  print(f'Num phospho sites:{len(phospho_sites):,} mapped non-primary IDs: {n_alt:,} missing IDs:{len(missing_ids):,}')
  
  #non_ref_dict = util.get_uniprot_columns(missing_ids, ['gene_primary', 'sequence'])
  
  return seq_dict
    
  # !phospho sites could be slightly difference in ref seqs
  

if __name__ == '__main__':
  
  csv_path = 'All_Total_phospho_Metacycle_RAIN_manual.csv'
 
  seq_dict = get_prot_sequences(csv_path)
  
  util.write_fasta('total_phospho_Mouse_seqs.fasta', seq_dict)
