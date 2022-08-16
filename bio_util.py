import re, os, multiprocessing, sys, time, gc
from subprocess import Popen, PIPE
from collections import defaultdict
import psutil

FASTA_SEQ_LINE = re.compile('(\S{59})(\S)')
FASTA_TRANS_TABLE = str.maketrans('JOBZU*', 'CPXXXX')
ALT_ID_URL = 'https://ftp.uniprot.org/pub/databases/uniprot/current%5Frelease/knowledgebase/complete/docs/sec%5Fac.txt'
UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/stream?query={}&format=tsv&fields={}'
HMMER_DIR = '/home/tjs23/hmmer-3.1b1-linux-intel-x86_64/binaries/'
NEWLINE_CHARS = 0
CLUSTAL_EXE = 'clustalo'
ALIGN_TMP_PATH = '__align_temp__'
PROG_CHARS = ('#', '-')
MAX_CORES = multiprocessing.cpu_count()
KEYBOARD_INTERRUPT_MSG = 'Parallel jobs stopped - KeyboardInterrupt'


# #  Job execution  # #    

def _parallel_job_wrapper(job, out_queue, target_func, data_item, args, kw):

  try:
    result = target_func(data_item, *args, **kw)
    out_queue.put((job, result), False)
 
  except KeyboardInterrupt as err: # Ignore in multiple sub-processes as this will be picked up only once later
    return
    

def parallel_run(target_func, job_data, common_args=(), common_kw={},
                 num_cpu=MAX_CORES, verbose=True, local_cpu_arg=None):
  # This does not use Pool.apply_async because of its inability to pass
  # pickled non-global functions
  
  from multiprocessing import Process, Manager
    
  job_data = [x for x in job_data if x]
  num_jobs = len(job_data)
  num_proc = min(num_cpu, num_jobs)
  procs = {} # Current processes
  queue = Manager().Queue() # Queue() # Collect output
  results = [None] * num_jobs
  gc.collect() # Mimimise mem footprint before fork
  
  if verbose:
    msg = 'Running {} for {:,} tasks on {:,} cores'
    info(msg.format(target_func.__name__, num_jobs, num_proc))
  
  k = 0
  prev = 0.0
  
  for j in range(num_jobs):
    
    if len(procs) == num_proc: # Full
      try:
        i, result = queue.get()
 
      except KeyboardInterrupt:
        critical(KEYBOARD_INTERRUPT_MSG)

      results[i] = result
      del procs[i]
      
      if verbose:
        k += 1
        f = k / float(num_jobs)
        
        if (f-prev) > 0.001: # Avoid very fine updates
          progress(k, num_jobs)
          prev = f
    
    if local_cpu_arg and (j >= num_proc): # After initial allocations
      arg, default = local_cpu_arg
      cpu_free = 1.0 - (psutil.cpu_percent()*0.01)
      kw_args = dict(common_kw)
      kw_args[arg] = max(default, int(cpu_free*num_proc*0.5))
    else:
      kw_args = common_kw
        
    args = (j, queue, target_func, job_data[j], common_args, kw_args)
    proc = Process(target=_parallel_job_wrapper, args=args)
    procs[j] = proc

    try:
      proc.start()
    
    except IOError:
      time.sleep(0.01)
      proc.start()

    except KeyboardInterrupt:
      critical(KEYBOARD_INTERRUPT_MSG)
      
    except Exception as err:
      raise(err)
      
  # Last waits
  
  while procs:
    try:
      i, result = queue.get()
      
    except KeyboardInterrupt:
      critical(KEYBOARD_INTERRUPT_MSG)
    
    results[i] = result
    del procs[i]
   
    if verbose:
      k += 1
      progress(k, num_jobs)
 
  if verbose:
    print('')
 
  #queue.close()
 
  return results


def get_num_cpu(min_cpu=2, max_frac=0.5, use_load=True):

  cpu_max = multiprocessing.cpu_count()
  
  if use_load:
    cpu_avail = 1.0 - (psutil.cpu_percent()*0.01)
  else:
    cpu_avail = 1.0
  
  num_cpu = max(min_cpu, int(cpu_avail*cpu_max*max_frac))
   
  return num_cpu
  

def report(msg, line_return=False):
  global NEWLINE_CHARS

  if line_return:
    fmt = '\r%%-%ds' % max(NEWLINE_CHARS, len(msg))
    sys.stdout.write(fmt % msg) # Must have enouch columns to cover previous msg
    sys.stdout.flush()
    NEWLINE_CHARS = len(msg)
    
  else: 
    if NEWLINE_CHARS:
      print('')
    print(msg)
    NEWLINE_CHARS = 0
    

def info(msg, line_return=False):

  report('INFO: ' + msg, line_return)
  

def progress(i, n):

  pc = (100.0*i)/n
  prog = int(pc)
  msg = '    |{}{}|{:7.2f}% [{:,}]'.format(PROG_CHARS[0] * prog, PROG_CHARS[1] * (100-prog), pc, n)
  report(msg, True)

 
def warn(msg):

  report('WARN: ' + msg)


def critical(msg):

  report('EXIT: ' + msg)
  report('STOP')
  sys.exit(0)


def get_uniprot_alt_ids(alt_id_url='https://ftp.uniprot.org/pub/databases/uniprot/current%5Frelease/knowledgebase/complete/docs/sec%5Fac.txt', cache_file='uniprot_sec_acc.txt'):

  from  urllib import request

  print(f'Fetching UniProt secondary accessions')
  
  alt_id_dict = {}
  
  if not os.path.exists(cache_file):
    req = request.Request(alt_id_url)
 
    with request.urlopen(req) as f:
       response = f.read().decode('utf-8')
       
       info(f' .. writing {cache_file}')
       with open(cache_file, 'w') as file_obj:
         file_obj.write(response)
  
  with open(cache_file) as file_obj:
    lines = file_obj.readlines()
 
  for i, line in enumerate(lines):
    if line.startswith('Secondary AC'):
      break
  
  for line in lines[i+1:]:
    if line.strip():
      sec_id, prim_id = line.split()
      alt_id_dict[sec_id] = prim_id
  
  return alt_id_dict
  
  
def _get_rand_string(size):
  
  import random, string
  
  return ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(size))
  
  
def _temp_fasta_file_path(named_seqs):

  file_path = '%s_%.3f%s%s' % (ALIGN_TMP_PATH, time.time(), _get_rand_string(8), '.fasta')
  write_fasta(file_path, named_seqs)
  
  return file_path


def fasta_string(named_seqs):

  if isinstance(named_seqs, dict):
    named_seqs = named_seqs.items()

  entries = [fasta_item(name, seq) for name, seq in named_seqs]

  return '\n'.join(entries)
  
  
def align_seqs_clustal(named_seqs, alignment=None, cpu_cores=None, aligned_input=False, full=False):
       
  if not cpu_cores:
    import multiprocessing
    cpu_cores = multiprocessing.cpu_count()

  if isinstance(named_seqs, dict):
    named_seqs = list(named_seqs.items())
  
  if alignment and  isinstance(alignment, dict):
    alignment = list(alignment.items())
  
  # Run Clustal Omega
  
  gen_cmd_args = ['--threads=%d' % cpu_cores]
  
  if full:
    gen_cmd_args += ['--full', '--full-iter', '--iter=2']
  
  if alignment:
    if not aligned_input and len(named_seqs) > 1:
      fasta_data = fasta_string(named_seqs).encode('utf-8')
      profile_temp_file = _temp_fasta_file_path(alignment)
      
      # Add multiple sequences to an existing alignment
       
      cmd_args = [CLUSTAL_EXE,
                 '--infile=-', # STDIN
                 '--profile1=%s' % profile_temp_file] + gen_cmd_args
                  
      proc = Popen(cmd_args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
      std_out_data, std_err_data = proc.communicate(fasta_data)
             
      if not std_out_data:
        # ClustalO doesn't like profiles that happen to have no gaps
        fasta_data = fasta_string(named_seqs + alignment)
        cmd_args = [CLUSTAL_EXE, '--infile=-'] + gen_cmd_args

        proc = Popen(cmd_args, stdin=PIPE, stdout=PIPE, stderr=PIPE)

        std_out_data, std_err_data = proc.communicate(fasta_data)
      
      os.unlink(profile_temp_file)
      
    else:
      profile_temp_file1 = _temp_fasta_file_path(named_seqs)
      profile_temp_file2 = _temp_fasta_file_path(alignment)
    
      # Add a single sequence to an existing alignment
      
      cmd_args = [CLUSTAL_EXE,
                 '--profile1=%s' % profile_temp_file1,
                 '--profile2=%s' % profile_temp_file2,
                 '--is-profile',  
                 '--threads=%d' % cpu_cores] # No further iteration (v. slow)

      proc = Popen(cmd_args, stdout=PIPE, stderr=PIPE)
      std_out_data, std_err_data = proc.communicate()

      if not std_out_data:
        # ClustalO doesn't like profiles that happen to have no gaps
        fasta_data = fasta_string(named_seqs + alignment).encode('utf-8')
        cmd_args = [CLUSTAL_EXE,
                   '--infile=-'] + gen_cmd_args
        
        proc = Popen(cmd_args, stdin=PIPE, stdout=PIPE, stderr=PIPE)

        std_out_data, std_err_data = proc.communicate(fasta_data)

      os.unlink(profile_temp_file1)
      os.unlink(profile_temp_file2)
  
  else:
    fasta_data = fasta_string(named_seqs).encode('utf-8')
    
    # Normal multiple alignment of input sequences
    
    cmd_args = [CLUSTAL_EXE,
               '--infile=-'] + gen_cmd_args
                 
    proc = Popen(cmd_args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    std_out_data, std_err_data = proc.communicate(fasta_data)
          
  std_out_data = std_out_data.decode('utf-8')
  std_err_data = std_err_data.decode('utf-8')
  
  # Read FASTA output

  aln_data = std_out_data.split('\n')
  align_dict = read_fasta(aln_data)
  
  if not align_dict:
    seq_names = ', '.join(['%s (%d)' % (x[0], len(x[1])) for x in named_seqs])
    msg1 = 'STDOUT: ' + std_out_data
    msg2 = 'STDERR: ' + std_err_data
    msg3 = 'Command used: {}'.format(' '.join(cmd_args))
    msg4 = 'Multiple alignment failed for sequences: {}'.format(seq_names)
    warn(msg1)
    warn(msg2)
    info(msg3)
    warn(msg4)
    return '\n'.join([msg1, msg2, msg3, msg4])
  
  return align_dict
  

def hmmer_search(query_seq_dict, database_fasta, iterations=1, 
                 cpu_cores=1, bit_score_cutoff=None, e_cutoff=None,
                 max_hits=None, matrix='BLOSUM62', query_hmm=None,
                 align_save_prefix=None, hmm_save_prefix=None, fast=False):
  
      
  query_fasta = '\n'.join([fasta_item(name, seq) for name, seq in query_seq_dict.items()]) + '\n'
  
  # Run jackhmmer (seq) or hmmsearch (HMM profile)
  
  if query_hmm:
    search_exe = 'hmmsearch'
  else:
    search_exe = 'jackhmmer'
    
  search_exe = os.path.join(HMMER_DIR, search_exe)
  
  cmd_args = [search_exe, '--cpu', str(cpu_cores), '-o', '/dev/null',
              '--domtblout', '/dev/stdout', '--noali']
  
  if fast:
    cmd_args += ['--F1', '0.01', '--F2', '0.001', '--F3', '1e-7'] # Improve speed
             
  if bit_score_cutoff:
    cmd_args += ['--domT', str(bit_score_cutoff)]
  
  elif e_cutoff:
    cmd_args += ['--domE', str(e_cutoff)]
  
  if query_hmm:
    cmd_args += [query_hmm, database_fasta]  
   
  else:
    cmd_args += ['--mx', matrix, '-N', str(iterations)]
    
    if align_save_prefix:
      cmd_args += ['--chkali', align_save_prefix]

    if hmm_save_prefix:
      cmd_args += ['--chkhmm', hmm_save_prefix]
      
    cmd_args += ['-', database_fasta] # seq STDIN
  
  try:
    proc = Popen(cmd_args, stdin=PIPE, stdout=PIPE)
 
  except Exception as err:
    warn('hmmer command failed')
    warn('Command used: "%s"' % ' '.join(cmd_args))
    warn(err)
    return
  
  std_out_data, std_err_data = proc.communicate(query_fasta.encode('utf-8'))
  
  if std_err_data:
    warn(cmd_args)
    warn(std_out_data, std_err_data)
    raise(Exception('HMMER fail'))
  
  lines = std_out_data.decode('utf-8').split('\n')
  results = defaultdict(list)
  
  for line in lines:
    if line and line[0] != '#':
      
      row = line.split()
      # name, acc, tlen, qname, qacc, qlen, evalue, bscore, bias, dom, ndom, ceval, ieval, dscore, dbias, hits-start, hit_end, q_start, q_end, acc, ...
      name = row[0]
      accession = row[1]
      
      if accession == '-':
        if '|' in name:
          accession = name.split('|')[1] # For UniProt
        else:
          accession = name
        
        name = ' '.join(row[22:])
      
      qacc = row[4]
      if qacc == '-':
        qacc = row[3]
        if '|' in qacc:
          qacc = qacc.split('|')[1] # For UniProt
      
      domain = int(row[9])
      bit_score = float(row[7])
      e_value = float(row[6])
      query_start = int(row[17])-1
      query_end = int(row[18])
      hit_start = int(row[15])-1
      hit_end = int(row[16])
      hit_len = int(row[2])
      
      result = [name, accession, e_value, bit_score, hit_len,
                query_start, query_end, hit_start, hit_end]
      
      if max_hits:
        if len(results[qacc]) < max_hits:
          results[qacc].append(result)    
      
      else:
        results[qacc].append(result)
        
          
  return results

def uniprot_header_acc(head):
  """ For  example
        tr|Q9WV19|Q9WV19_MOUSE Cytochrome P450 OS=Mus musculus OX=10090 GN=Cyp2g1 PE=2 SV=1 
      becomes
        Q9WV19
  """
  
  return head.split('|')[1]
  

def read_fasta(stream_or_path, as_dict=True, head_processor=None):
  
  if isinstance(stream_or_path, str):
    stream = open(stream_or_path)
  elif isinstance(stream_or_path, (list, tuple)):
    stream = stream_or_path
  else:
    stream = stream_or_path
  
  named_seqs = []
  append = named_seqs.append
  name = None
  seq = []

  for line in stream:
    line = line.strip()
    
    if not line:
      continue
    
    if line[0] == '>':
      if name:
        append((name, ''.join(seq).translate(FASTA_TRANS_TABLE)))

      seq  = []
      name = line[1:]
      
      if head_processor:
        name = head_processor(name)
    
    elif name:
      seq.append(line)

  if name:
    append((name, ''.join(seq).translate(FASTA_TRANS_TABLE)))

  if as_dict:
    return dict(named_seqs)
  else:
    return named_seqs


def get_uniprot_columns(queries, columns=['id','protein_name','gene_primary'], batch_size=100):
  
  from  urllib import request
  
  #info(f'Fetching UniProt data for {len(queries):,} entries')
  
  a = 0
  n = len(queries)
  
  uniprot_dict = {}
  
  while a < n:
    b = min(n, a+batch_size)
    info(f' .. {a:,} - {b:,}', line_return=True)
    
    query_ids= '%20OR%20'.join(list(queries[a:b]))
    fields = '%2C'.join(['accession'] + columns)
    req = request.Request(UNIPROT_URL.format(query_ids, fields))
 
    with request.urlopen(req) as f:
       response = f.read()
 
    lines = response.decode('utf-8').split('\n')
    
    for line in lines[1:]:
      line = line.rstrip('\n')
 
      if not line:
        continue
        
      row = line.split('\t')  
      uniprot_dict[row[0]] = row[1:]
    
    a = b 
     
  info(f' .. done')
  
  if len(columns) == 1:
    for pid in uniprot_dict:
      uniprot_dict[pid] = uniprot_dict[pid][0]
  
  return uniprot_dict


def fasta_item(name, seq, end=''):

  seq = seq.replace(u'\ufffd', '')
  seq = seq.upper()
  seq = re.sub('\s+','',seq)
  seq = seq[0] + FASTA_SEQ_LINE.sub(r'\1\n\2',seq[1:])
 
  return '>%s\n%s%s' % (name, seq, end) 


def write_fasta(file_path, named_seqs):
  
  if isinstance(named_seqs, dict):
    named_seqs = list(named_seqs.items())
  
  n = 0
  with open(file_path, 'w') as file_obj:
    write = file_obj.write

    for name, seq in named_seqs:
      write('%s\n' % fasta_item(name, seq) )
      n += 1
      
  info(f'Written {n:,} sequences to {file_path}')
  
