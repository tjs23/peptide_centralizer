import bio_util as util
from glob import glob
from collections import defaultdict

def reciprocal_filter(homologue_path, orthologue_path, ref_proteome_fasta):
  print(f'Filtering reciprocal hits in {homologue_path}')
  orthologue_dict = {}
  
  ref_seqs = util.read_fasta(ref_proteome_fasta, head_processor=util.uniprot_header_acc)
  
  with open(homologue_path) as file_obj:
    for line in file_obj:
      query, *hids = line.rstrip().split()
      oid, qid = query.split(':')
      
      seq_dict = util.get_uniprot_columns(hids, ['sequence'])
      
      print(f' .. checking {oid} {qid} : {len(hids)} homologues')

      hit_dict = util.hmmer_search(seq_dict, ref_proteome_fasta, iterations=1, cpu_cores=10, bit_score_cutoff=None,
                                   e_cutoff=1e-20, max_hits=3, matrix='BLOSUM62', query_hmm=None,
                                   align_save_prefix=None, hmm_save_prefix=None, fast=False) or []
      
      rid = None
      
      for pid in seq_dict:
        if pid in ref_seqs:
          rid = pid
          print(f' .. using reciprocation ref sequence {rid}')
          break  
      
      if not rid:
       print(f' .. reciprocation ref missing, reverting to {qid}')
       rid = qid
      
      oids = [rid]
      
      for pid in hit_dict:
        if pid == rid:
          continue
      
        hits = hit_dict[pid]
        
        if hits:
          best_score = hits[0][3]
          
          for name, accession, e_value, bit_score, hit_len, query_start, query_end, hit_start, hit_end in hits:
            if bit_score < 0.999 * best_score:
              break
            
            if accession in (rid, qid):
              oids.append(pid)
              break
      
      print(f' .. checking {qid}:{rid} : confirmed {len(oids)}')
      orthologue_dict[query] = oids

  with open(orthologue_path, 'w') as file_obj:
    for query in orthologue_dict:
      oids = ' '.join(sorted(orthologue_dict[query]))

      file_obj.write(f'{query} {oids}\n')
  
  print(f'Written {orthologue_path}')
  

def gather_homologues(query_path, proteome_fastas, homologue_path):
  
  seq_dict = util.read_fasta(query_path)
  homologues = defaultdict(list)
  
  n = len(seq_dict)
  
  for database_fasta in proteome_fastas:
    util.info(database_fasta)
    
    hit_dict  = util.hmmer_search(seq_dict, database_fasta, iterations=1, cpu_cores=10, bit_score_cutoff=None,
                                  e_cutoff=1e-20, max_hits=1, matrix='BLOSUM62', query_hmm=None,
                                  align_save_prefix=None, hmm_save_prefix=None, fast=False) or []
          
    for query in hit_dict:
      hits = hit_dict[query]
   
      for name, accession, e_value, bit_score, hit_len, query_start, query_end, hit_start, hit_end in hits:
        homologues[query].append(accession)           
        break    

  with open(homologue_path, 'w') as file_obj:
    
    for query in homologues:
      oids = ' '.join(sorted(homologues[query]))

      file_obj.write(f'{query} {oids}\n')
      
  print(f'Written {homologue_path}')
  

if __name__ == '__main__':
  
  query_path = 'total_phospho_Mouse_seqs.fasta'
  
  ref_proteome_fasta = 'proteomes/UP000000589_Mus_musculus.fasta'
  
  proteome_fastas = glob('proteomes/UP*.fasta')
  
  homologue_path = 'total_phospho_Mouse_homologues.txt'

  orthologue_path = 'total_phospho_Mouse_orthologues.txt'
  
  gather_homologues(query_path, proteome_fastas, homologue_path)
  
  reciprocal_filter(homologue_path, orthologue_path, ref_proteome_fasta)
