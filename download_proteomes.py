from  urllib import request

PROTEOME_URL = 'https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco&format=tsv&query=(taxonomy_name%3A{})%20AND%20(proteome_type%3A1)'
FASTA_URL = 'https://rest.uniprot.org/uniprotkb/stream?query=(proteome:{})&format=fasta'

def download_uniprot_proteome_clade_fastas(tax_clade, url=PROTEOME_URL, top=None):
    
    furl = url.format(tax_clade)

    print(f'Proteome query {furl}')

    with request.urlopen(request.Request(furl)) as f:
       lines = f.read().decode('utf-8').split('\n')
    
    proteomes = []
    
    for line in lines[1:]:
      line = line.rstrip('\n')
      
      if not line:
        continue
      
      upid, species, tax_id, count, busco = line.split('\t')
      
      completeness = float(busco[busco.index('C:')+2:busco.index('%[')])
      proteomes.append([completeness, upid, species, tax_id])
    
    proteomes.sort(reverse=True)
    
    if top:
      proteomes = proteomes[:top]
   
    for completeness, upid, species, tax_id in proteomes:
      print(f' .. {completeness}% {upid} {species} {tax_id}')
      
    return proteomes
  
  
def download_uniprot_proteome_fasta(species, uid, url=FASTA_URL, overwrite=False):
    
    furl = url.format(uid)
    
    print(f'Query {species}')
    
    print(f' .. download {furl}')
       
    out_file_path = f'{uid}_{species}.fasta'
     
    req = request.Request(furl)
 
    with request.urlopen(req) as f:
       response = f.read()
 
    lines = response.decode('utf-8')

    print(f' .. obtained {lines.count(">"):,} sequences')
    

    with open(out_file_path, 'w') as file_obj:
      file_obj.write(lines)

    print(f' .. saved {out_file_path} ')


if __name__ == '__main__':
  
  proteomes = []
  proteomes += download_uniprot_proteome_clade_fastas('Primates', top=12)
  proteomes += download_uniprot_proteome_clade_fastas('Scandentia')
  proteomes += download_uniprot_proteome_clade_fastas('Dermoptera')
  proteomes += download_uniprot_proteome_clade_fastas('Glires')
  
  for completeness, upid, species, tax_id in proteomes:
    species = species.split(' (')[0].replace(' ','_')
    download_uniprot_proteome_fasta(species, upid)
