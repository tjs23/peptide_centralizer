import bio_util as util


def align_orthologues(query_fasta_path, orthologue_list_path, alignments_path):
  
  qseq_dict = util.read_fasta(query_fasta_path)
  
  print(f'Aligning sequences with IDs from {orthologue_list_path}')
  
  align_inputs = []
  singular_inputs = [] 
  
  i = 0
  orthologues = []
  all_pids = set()
  with open(orthologue_list_path) as file_obj:
    for line in file_obj:
      query, *pids = line.split()
      orthologues.append((query, pids))
      all_pids.update(pids)
  
  all_pids = sorted(all_pids)
  all_seq_dict = util.get_uniprot_columns(all_pids, ['sequence'])
  queries = []
  
  for query, pids in orthologues:
    seq_dict = {x:all_seq_dict[x] for x in pids}
    qid, pid = query.split(':')
  
    if pid not in seq_dict: # Primary ID is not a reference sequence; no big deal
      seq_dict[pid] = qseq_dict[query]
      
    if len(seq_dict) < 2:
      singular_inputs.append((query, seq_dict))
    else:
      align_inputs.append(seq_dict)
      queries.append(query)
    
    util.info(f' .. {i:,}', line_return=True)
    i += 1
      
  print('')
  
  align_dicts = util.parallel_run(util.align_seqs_clustal, align_inputs, common_kw={'cpu_cores':2})
  align_dicts = list(zip(queries, align_dicts))
  align_dicts += singular_inputs
 
  with open(alignments_path, 'w') as out_file_obj:
    write = out_file_obj.write
        
    i = 0
    for query, align_dict in align_dicts:
      fasta_data = []
      for aid, aseq in align_dict.items():
        head = f'{query}|{aid}'
        fasta_data.append(util.fasta_item(head, aseq))
 
      write('\n'.join(fasta_data) + '\n')
 
      i += 1
      util.info(f' .. {i:,}', line_return=True)
 
    util.info(f' Written {i:,} alignments to {alignments_path}')
 
  
      
if __name__ == '__main__':
  
  query_fasta_path = 'total_phospho_Mouse_seqs.fasta'
  
  orthologue_list_path = 'total_phospho_Mouse_orthologues.txt'
  
  alignments_path = 'total_phospho_Mouse_alignments.fasta'
  
  align_orthologues(query_fasta_path, orthologue_list_path, alignments_path)
