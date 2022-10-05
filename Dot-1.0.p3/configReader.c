#include "configReader.h"

struct config *command_line_parser (int argc, char* argv[]) {
  struct config *cfg;
  int c;
  int opt_index = 0;
  static struct option long_options[] = {
    {"help",      no_argument,       0,  'h' },
    {"config",    required_argument, 0,  'c' },
    {"sequence",    required_argument, 0,  's' },
    {"minmotif",    required_argument, 0,  'l' },
    {"maxmotif",    required_argument, 0,  'L' },
    {"minmatch",    required_argument, 0,  'm' },
    {"output",    required_argument, 0,  'o' },
    {"maxinsert",    required_argument, 0,  'I' },
    {"maxgaps",   required_argument, 0,  'G' },
    {"threads",   required_argument, 0,  't' },
    {"verbose",   no_argument, 0,  'v' },
    {"version",   no_argument, 0,  'V' },
    {0,           0,                 0,  0   }
  };
  cfg = (struct config *) malloc (sizeof (struct config));
  if (cfg == NULL) {
    perror ("Memory allocation failure\n");
    return NULL;
  }
  cfg->flags = 0;
  cfg->thread_max = 1;
  cfg->filter_type = NONE;
  cfg->allow_overlap = true;
  cfg->fvalue = 1;  /*  MinMatch   */
  cfg->nvalue = 2;  /*  min motif length  */
  cfg->xvalue = 30; /*  max motif length  */
  cfg->jvalue = 0;   /*  MaxInsert - XXX called max jumps */
  cfg->gvalue = 0;   /*  MaxGaps  */
  cfg->tollerance = 0;
  cfg->min_TR_len = 12;
  cfg->min_purity = 0.1;
  cfg->verbose = false;
  /*  Parameters from command line  */
  while ((c = getopt_long(argc, argv,"hc:s:l:L:m:o:I:G:t:vV", 
			  long_options, &opt_index )) != -1) {
    switch (c) {
    case 'v':  /*  enable verbose output  */
      cfg->verbose = true;
      break;
    case 'c':  /*  config file name  */
      cfg->flags += C_FLAG;
      strcpy (cfg->cvalue, optarg);
      break;
    case 's':  /*  input seqyuence file name  */
      cfg->flags += S_FLAG;
      strcpy (cfg->svalue, optarg);
      break;
    case 'l':   /*  min motif length  */
      cfg->flags += N_FLAG;
      cfg->nvalue = atoi (optarg);
      break;
    case 'L':   /*  max motif length  */
      cfg->flags += X_FLAG;
      cfg->xvalue = atoi (optarg);
      break;
    case 'm':   /*  MinMatch  */
      cfg->flags += F_FLAG;
      cfg->fvalue = (MATCH_ARRAY_TYPE) atof (optarg);
      break;
    case 'o':   /*  Output file name  */
      cfg->flags += O_FLAG;
      strcpy (cfg->output_filename, optarg);      
      break;
    case 'h':  /*  print help and exit  */
      cfg->flags += U_FLAG;
      break;
    case 'I':   /*  MaxInsert  */
      cfg->flags += J_FLAG;
      cfg->jvalue = atoi (optarg);
      break;
    case 'G':  /*  MaxGaps  */
      cfg->flags += G_FLAG;
      cfg->gvalue = atoi (optarg);
      break;
    case 't':    /*  Number of threads  */
      cfg->flags += T_FLAG;
      cfg->thread_max = atoi (optarg);
      break;
    case 'V':
      print_version ();
      free (cfg);
      exit (EXIT_SUCCESS);
    case '?':
      print_usage ();
      free (cfg);
      exit (EXIT_FAILURE);
    default:
      free (cfg);
      exit (EXIT_FAILURE);
    }
  }
  return cfg;
}


MATCH_ARRAY_TYPE **read_weights_matrix (struct paramList *param) {
  struct paramList *elem;
  MATCH_ARRAY_TYPE **m;
  int i;
  m = (MATCH_ARRAY_TYPE **) malloc (sizeof (MATCH_ARRAY_TYPE *) * MAXIND_ARRAY);
  if (m == NULL) {
    perror ("Memory allocation failure reading weights\n");
    return NULL;
  }
  m[0] = (MATCH_ARRAY_TYPE *) malloc (sizeof (MATCH_ARRAY_TYPE) * MAXIND_ARRAY * MAXIND_ARRAY);
  if (m[0] == NULL) {
    perror ("Memory allocation failure reading weights\n");
    return NULL;
  }
  /*  Set pointers  */
  for (i = 0; i < MAXIND_ARRAY; i ++) m[i] = m[0] + (MAXIND_ARRAY * i);
  /*  Parameters from config file */
  for (elem = param; elem != NULL; elem = elem->next) {
    if (strlen (elem->param) != 2) continue;
    if (elem->param[0] != 'A' && elem->param[0] != 'C' && \
	elem->param[0] != 'G' && elem->param[0] != 'T' &&   \
	elem->param[0] != 'N' && elem->param[0] != 'a' &&   \
	elem->param[0] != 'c' && elem->param[0] != 'g' &&   \
	elem->param[0] != 't' && elem->param[0] != 'n' ) continue;
    if (elem->param[1] != 'A' && elem->param[1] != 'C' && \
	elem->param[1] != 'G' && elem->param[1] != 'T' &&   \
	elem->param[1] != 'N' && elem->param[1] != 'a' &&   \
	elem->param[1] != 'c' && elem->param[1] != 'g' &&   \
	elem->param[1] != 't' && elem->param[1] != 'n' ) continue;

    m[(int) elem->param[0]][(int) elem->param[1]] = atof (elem->value);
    m[(int) elem->param[1]][(int) elem->param[0]] = m[(int) elem->param[0]][(int) elem->param[1]]; 
  }
  return m;
}

void free_weights_matrix (MATCH_ARRAY_TYPE **m) {
  free (m[0]);
  free (m);
  m = NULL;
}



void param_list_parser (struct paramList *par_list, struct config *cfg) {
  struct paramList *elem;
  int i;
  /** Parameters from config file */
  for (elem = par_list; elem != NULL; elem = elem->next) {
    if (strlen (elem->param) == 2) continue;
    /*  Convert to lowercase to be less fussy with param detection  */
    for(i = 0; elem->param[i] != '\0'; i ++) {
      elem->param[i] = tolower(elem->param[i]);
    }
	    
    if (strcmp (elem->param, "filtertype") == 0) {
      if ((cfg->flags & FILT_FLAG) == 0) {
	cfg->flags += FILT_FLAG;
	/*  Convert to lowercase to be less fussy with param detection  */
	for(i = 0; elem->value[i] != '\0'; i ++) {
	  elem->value[i] = tolower(elem->value[i]);
	}      
	if (strcmp (elem->value, "heavy") == 0) cfg->filter_type = HEAVY;
	if (strcmp (elem->value, "fair") == 0) cfg->filter_type = FAIR;
	if (strcmp (elem->value, "light") == 0) cfg->filter_type = LIGHT;
	if (strcmp (elem->value, "threshold") == 0) cfg->filter_type = THRESHOLD;
	if (strcmp (elem->value, "none") == 0) cfg->filter_type = NONE;
      }
    }
    if (strcmp (elem->param, "tolerance") == 0) {  
      if ((cfg->flags & TOLL_FLAG) == 0) {
	cfg->flags += TOLL_FLAG;
	cfg->tollerance = (float) atof (elem->value);
      }
    }
    if (strcmp (elem->param, "allowoverlap") == 0) {
      if ((cfg->flags & OVER_FLAG) == 0) {
	cfg->flags += OVER_FLAG;
	if (elem->value[0] == 'Y') cfg->allow_overlap = true;
	if (elem->value[0] == 'N') cfg->allow_overlap = false;
	if (elem->value[0] == 'y') cfg->allow_overlap = true;
	if (elem->value[0] == 'n') cfg->allow_overlap = false;
      }
    }
    if (strcmp (elem->param, "mintrlen") == 0) {
      if ((cfg->flags & MTRL_FLAG) == 0) {
	cfg->flags += MTRL_FLAG;
	cfg->min_TR_len = atoi (elem->value);
      }
    }
    if (strcmp (elem->param, "minpurity") == 0) {
      if ((cfg->flags & MPUR_FLAG) == 0) {
	cfg->flags += MPUR_FLAG;
	cfg->min_purity = (float) atof (elem->value);	
      }
    }
    
    if (strcmp (elem->param, "minmotiflen") == 0) {
      if ((cfg->flags & N_FLAG) == 0) {
	cfg->flags += N_FLAG;
	cfg->nvalue = atoi (elem->value);
      }
    }
    if (strcmp (elem->param, "maxmotiflen") == 0) {
      if ((cfg->flags & X_FLAG) == 0) {
	cfg->flags += X_FLAG;
	cfg->xvalue = atoi (elem->value);
      }
    }
    if (strcmp (elem->param, "minmatch") == 0) {
      if ((cfg->flags & F_FLAG) == 0) {
	cfg->flags += F_FLAG;
	cfg->fvalue = (MATCH_ARRAY_TYPE) atof (elem->value);
      }
    }
    if (strcmp (elem->param, "sequence") == 0) {
      if ((cfg->flags & S_FLAG) == 0) {
	cfg->flags += S_FLAG;
	strcpy (cfg->svalue, elem->value); 
      }
    }
    if (strcmp (elem->param, "outfile") == 0) {
      if ((cfg->flags & O_FLAG) == 0) {
	cfg->flags += O_FLAG;
	strcpy (cfg->output_filename, elem->value); 
      }
    }
    if (strcmp (elem->param, "maxinsert") == 0) {
      if ((cfg->flags & J_FLAG) == 0) {
	cfg->flags += J_FLAG;
	cfg->jvalue = atoi (elem->value);
      }
    }
    if (strcmp (elem->param, "maxgaps") == 0) {
      if ((cfg->flags & G_FLAG) == 0) {
	cfg->flags += G_FLAG;
	cfg->gvalue = atoi (elem->value);
      }
    }
    if (strcmp (elem->param, "treads") == 0) {
      if ((cfg->flags & T_FLAG) == 0) {
	cfg->flags += T_FLAG;
	cfg->thread_max = atoi (elem->value);
      }
    }
  }
}


void free_paramList (struct paramList* pl) {
  struct paramList *elem_succ = pl, *curr = NULL;
  while (elem_succ != NULL) {
    curr = elem_succ;
    elem_succ = elem_succ->next;
    free (curr);
  }
  pl = NULL;	
}

struct paramList *loadConfigFromFile (char *filename) {
  struct paramList *start = NULL, *l = NULL;
  char *key, *val;
  FILE *pf;
  char line [MAX_FILELINE]; 
  int i, lineNumber = 0;
  if ((pf = fopen (filename, "r" )) == NULL) {
    printf ("Unable to open config file\n");
    return NULL;
  }
  while (fgets (line, sizeof line, pf) != NULL ) {
    lineNumber ++;
    /*  remove padding at the end of the line  */
    for (i = strlen (line) - 1; i >= 0; i --) {
      if (line[i] == ' ' || line[i] == '\n' || line[i] == '\t' || line[i] == '\r') line[i] = '\0';
      else break;
    }
    for (i = 0; line[i] != '\0'; i ++) {  /*  skip white spaces  */
      if (line[i] != ' ') break;
    }
    if (line[i] == '\0') continue;  /*  it is an empty line  */
    if (line[i] == '#') continue;   /*  it is a comment   */
    val = line + i;  /*   set the offset of the beginning of the string remobing spaces  */
    key = strsep (&val, "=");
    if (val == NULL) {
      printf ("Error in line %d: this is not a key = value pair\n", lineNumber);
      return NULL;
    }
    /*  remove trailing spaces from key and val  */
    while (val[0] != '\0') {
      if (val[0] == ' ') val ++;
      else break;
    }
    if (val[0] == '\0') {
      printf ("Error in line %d: missing value\n", lineNumber);
      return NULL;
    }
    for (i = strlen (key) - 1; i >= 0; i --) {
      if (key[i] == ' ' || key[i] == '\t') key[i] = '\0';
      else break;
    }
    if (key[0] == '\0') {
      printf ("Error in line %d: missing key\n", lineNumber);
      return NULL;
    }
    if (start == NULL) {
      start = (struct paramList *) malloc (sizeof (struct paramList));
      l = start;
    }
    else {
      l->next = (struct paramList *) malloc (sizeof (struct paramList));
      l = l->next;
    }
    if (l == NULL) return NULL;  /*  Does not free memory - the program will terminate  */
    strncpy (l->param, key, MAX_FILELINE);
    strncpy (l->value, val, MAX_FILELINE);
    /*  printf ("*%s*, *%s*\n",key, val);   */
  }
  fclose (pf);
  return start;
}
