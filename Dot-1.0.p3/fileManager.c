#include "fileManager.h"

void __filemanager_read_id (char *filename, char *id_buff) {
  int i, j, fn_len, start_pos, end_pos;
  fn_len = strlen (filename);
  start_pos = 0;
  /*  Skip the possiblw path in the filename  */
  for (i = 0; i < fn_len; i ++) {
    if (filename[i] == '/' || filename[i] == '\'') {  /*  Both windows and linux  */
      start_pos = i + 1;
    }
  }
  /* skip the file extension  */
  for (end_pos = fn_len - 1; end_pos > start_pos && filename[end_pos] != '.'; end_pos --);
  /*  Copy file name bounded to the maximum label length  */
  for (i = start_pos, j = 0; i < end_pos && j < MAX_LABEL_LENGTH - 1; i ++, j ++) {
    id_buff[j] = filename[i];
  }
  id_buff[j] = '\0';
}

seqfile_t __filemanager_get_filetype (struct filemanager *fmobj) {
  unsigned int i;
  for (i = fmobj->offset; i < fmobj->buffer_size; i ++) {
    switch (fmobj->buffer[i]) {
    case '>':
      return FASTA;
    case '@':
      return FASTQ;
    case ' ':
      break;
    default:
      return UNKNOWN;
    }
  }
  return UNKNOWN;
}
void filemanager_destroy (struct filemanager *fmobj) {
  fclose ( fmobj->pf );
  free (fmobj);
  fmobj = NULL;
}

struct filemanager *filemanager_init (char *filename) {
  struct filemanager *fmobj;
  if ((fmobj = (struct filemanager *) malloc (sizeof (struct filemanager))) == NULL) {
    perror ("Memory allocation failure parsing input\n");
    return NULL;
  }
  if ((fmobj->pf = fopen (filename, "r")) == NULL) {
    free (fmobj);
    perror ("Error opening input file\n");
    return NULL;
  }
  __filemanager_read_id (filename, fmobj->empty_identifier);
  fmobj->offset = 0;
  fmobj->buffer_size = fread (fmobj->buffer, 1, BUFF_SIZE, fmobj->pf);
  /*  Check reading error  */
  if (fmobj->buffer_size < BUFF_SIZE && feof (fmobj->pf) == 0) {
    perror ("Error reading input file\n");
    free (fmobj);
    return NULL;
  }
  fmobj->filetype = __filemanager_get_filetype (fmobj);
  if (fmobj->filetype == UNKNOWN) {
    perror ("Input file format not recognized\n");
    free (fmobj);
    return NULL;
  }
  return fmobj;
}

struct sequence_t *__filemanager_next_seq (struct filemanager *fmobj, struct sequence_t *seq) {
  char next_char;
  parse_status_t parse_status = H_PRE_SI;
  int qual_read_char = 0;  /*  Count of the number of caracters read in quality score of fastq */
  seq->label_size = 0;
  seq->sequence_size = 0;
  if (fmobj->finish == true) {
    free (seq->sequence);
    free (seq);
    seq = NULL;
    return NULL;
  }
  while (1) {
    if (fmobj->offset == fmobj->buffer_size) {
      if (feof (fmobj->pf) == 1) {  /*  End of file and end of buffer  */
	fmobj->finish = true;
	break;
      }
      fmobj->offset = 0;
      fmobj->buffer_size = fread (fmobj->buffer, 1, BUFF_SIZE, fmobj->pf);
      /*  Check reading error  */
      if (fmobj->buffer_size < BUFF_SIZE && feof (fmobj->pf) == 0) {
	perror ("Error reading input file\n");
	free (seq->sequence);
	free (seq);
	seq = NULL;
	return NULL;
      }
    }
    next_char = fmobj->buffer[fmobj->offset];
    fmobj->offset ++;
    /*  From here on I process the sequence/header  */
    switch (parse_status) {
    case H_PRE_SI:
      switch (next_char) {
      case '>':
      case '@':
	parse_status = H_PRE_LABEL;
	break;
      case ' ':
	break;
      default:
	perror ("Input file parsing error");
	free (seq->sequence);
	free (seq);
	seq = NULL;
	return NULL;
      }
      break;
    case H_PRE_LABEL:
      switch (next_char) {
      case '\n':
	strcpy (seq->label, fmobj->empty_identifier);
	seq->label_size = strlen (fmobj->empty_identifier);
	parse_status = SEQUENCE;
	break;
      case ' ':
	break;
      default:
	parse_status = H_LABEL;
	fmobj->offset --;  /*  Reprocess this caracter  */
	break;
      }
      break;
    case H_LABEL:
      switch (next_char) {
      case '\n':
	parse_status = SEQUENCE;
	break;
      case ' ':
	/*  without break to allow spaces in fastq (instead trim in fasta)  */
	if (fmobj->filetype == FASTA) parse_status = H_POST_LABEL;
      default:
	if (seq->label_size < MAX_LABEL_LENGTH - 1) {  /*  Prevent exceed buffer size with \0   */
	  seq->label[seq->label_size] = next_char;
	  seq->label_size ++;
	  seq->label[seq->label_size] = '\0';
	}
	break;
      }
      break;
    case H_POST_LABEL:
      if (next_char == '\n') parse_status = SEQUENCE;
      break;
    case SEQUENCE:
      switch (next_char) {
      case ' ':
	break;
      case '\n':
	break;
      case '>':  /*  No break because in fastq the caracter is evaluated  */
	if (fmobj->filetype == FASTA) {
	  fmobj->offset --;  /*  Reprocess this caracter  */
	  return seq;
	}
      case '+':  /*  No break because in fasta the caracter is evaluated  */
	if (fmobj->filetype == FASTQ) {
	  parse_status = FQ_PLUS;
	  break;
	}
      default:
	/*   Buffer extension is required here  */
	if (seq->sequence_size == seq->buffer_size - 1) {  
	  seq->sequence = (char *) realloc (seq->sequence, (seq->buffer_size + BUFF_SIZE) * sizeof (char));
	  if (seq->sequence == NULL) {
	    perror ("Error reallocating memory while reading the input\n");
	    free (seq);
	    seq = NULL;
	    return NULL;
	  }
	  seq->buffer_size +=  BUFF_SIZE;
	}
	seq->sequence[seq->sequence_size] = next_char;
	seq->sequence_size ++;
	seq->sequence[seq->sequence_size] = '\0';
	break;
      }
      break;
    case FQ_PLUS:
      if (next_char == '\n') parse_status = FQ_SCORE;
      break;
    case FQ_SCORE:
      switch (next_char) {
      case '\n':
	break;
      case '@':  /*  without break - if it is not the new element it is part of the score  */
	if (qual_read_char >= seq->sequence_size) {
	  fmobj->offset --;  /*  Reprocess this caracter  */
	  return seq;
	}
      default:
	qual_read_char ++;
	break;
      }
      break;
    }
  }
  return seq;
}

struct sequence_t *filemanager_next_seq (struct filemanager *fmobj, struct sequence_t *seq) {
  if (seq == NULL) {
    seq = (struct sequence_t *) malloc (sizeof (struct sequence_t));
    if (seq == NULL) {
      perror ("Memory error reading a new sequence\n");
      return NULL;
    }
    seq->sequence = (char *) malloc (BUFF_SIZE * sizeof (char));
    if (seq->sequence == NULL) {
      free (seq);
      perror ("Memory error reading a new sequence\n");
      return NULL;
    }
    seq->buffer_size = BUFF_SIZE;
  }
  seq->label_size = 0;
  seq->sequence_size = 0;
  if (fmobj->filetype == FASTQ || fmobj->filetype == FASTA) 
    return __filemanager_next_seq (fmobj, seq);
  return NULL;
}
