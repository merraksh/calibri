/*
 * Name:    cmdline.cpp
 * Author:  Pietro Belotti
 * Purpose: procedures for reading/managing command line arguments
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cmdLine.hpp"

#define MALLOC_BLOCK 20


#ifdef _MSC_VER
#define strcpy strcpy_s
#define sprintf sprintf_s
#else
#define strcpy(a,b,c) strncpy(a,c,b)
#endif

/*
 * Action to be taken for each parameter
 */

int act (char ***argv, int *argc, enum etype tact, void *par) {

  switch (tact) {

  case TTOGGLE:                                       if (par) *(char   *) par = 1;             break; /**< Toggle some parameter */
  case TCHAR:   if (!((*argc)--)) return -1; ++*argv; if (par) *(char   *) par =      ***argv;  break; /**< Char param */
  case TINT:    if (!((*argc)--)) return -1; ++*argv; if (par) *(int    *) par = atoi (**argv); break; /**< Integer param */
  case TDOUBLE: if (!((*argc)--)) return -1; ++*argv; if (par) *(double *) par = atof (**argv); break; /**< Double param */
  case TSTRING: if (!((*argc)--)) return -1; ++*argv;

    if (par) {
      *(char **) par = (char *) realloc (*(char **) par, (1 + strlen (**argv)) * sizeof (char));
      strcpy (* (char **) par, 1 + strlen (**argv), **argv);  /**< String param */
      //printf ("read param %s\n", *par);
    }

    break;

  default: printf ("switch{} mismatch in act ()\n"); return -1;
  }

  return 0;
}


/*
 * Set default value of the command line arguments
 */

void set_default_args (tpar *options) {

  for (; options -> shortopt; options++) 

    if (options -> par)

      switch (options -> type) {

      case TINT:    * (int    *) (options -> par) = (int)    options -> initial; break;
      case TTOGGLE: * (char   *) (options -> par) = (char)   options -> initial; break;
      case TCHAR:   * (char   *) (options -> par) = (char)   options -> initial; break;
      case TDOUBLE: * (double *) (options -> par) = (double) options -> initial; break;
      case TSTRING: strcpy (*(char **) options -> par, 1, "");                      break;
      default: ;
      }
}


/*
 * Read argument line
 */ 

char **readargs (int argc, char **argv, tpar *options) {

  char **filenames = NULL;

  char *cumul;
  char *pname = *argv;

  int nfiles = 0;
  tpar *p2;

  /*
   * parse arguments, set corresponding parameters
   */

  while (--argc) {

    ++argv;

    if ((options) && 
	(** argv == '-') &&    /* arg starts with a "-"                               */
	((*argv) [1])) {       /* not a single "-" (meaning I should read from stdin) */

      cumul = *argv + 1;

      // - long option ----------------------------------------------

      if (*cumul == '-') {  /** long option  ("--something [val]") */

	++cumul;

	for (p2 = options; p2 -> shortopt; p2++) {

	  if (!strcmp (cumul, p2 -> longopt)) {

	    //if (p2 -> type != TTOGGLE) 
	    //  ++argv;		  

	    if (//p2 -> par &&
		act (&argv, &argc, p2->type, (p2->par))) {

	      printf ("%s: missing value", *(argv-1)); 
	      return NULL;
	    }

	    break;
	  }
	}

	if (!(p2 -> shortopt)) {
	  printf ("%s: unrecognized option \"--%s\"\n", pname, cumul); 
	}
      }

      // - short option ---------------------------------------------

      else {                   /** short option ("-s [val]") */

	for (; *cumul; cumul++) {

	  for (p2 = options; p2 -> shortopt; p2++)

	    if (*cumul == p2->shortopt) { 

	      if (act (&argv, &argc, p2->type, (p2->par))) {
		printf ("%s: missing value", *(argv-1)); 
		return NULL;
	      }
	      break;
	    }

	  if (!(p2 -> shortopt)) {

	    printf ("%s: unrecognized option \"-%c\"\n", pname, *cumul); 
	  }
	}
      }
    }
    else {

      if (!(nfiles % MALLOC_BLOCK)) 
	filenames = (char **) realloc (filenames, (nfiles + MALLOC_BLOCK) * sizeof (char *));

      filenames [nfiles] = (char *) malloc ((strlen (*argv) + 1) * sizeof (char));

      strcpy (filenames [nfiles++], 1 + strlen (*argv), *argv);
    }
  }

  if (filenames) {

    if (!(nfiles % MALLOC_BLOCK)) 
      filenames = (char **) realloc (filenames, (nfiles + 1) * sizeof (char *));

    filenames [nfiles] = NULL;
  }

  return filenames;
}


/*
 * Print help message
 */

void print_help (char *pname, tpar *opt) {

  int i;

  char arg [5];
  char helpline [80];

  printf ("Usage: %s [options] file...\nOptions are:\n", pname);

  for (i=0; opt [i]. shortopt; i++) {

    switch (opt [i]. type) {
    case TTOGGLE: *arg=0;                break; /*sprintf (arg, "    "); break;*/
    case TCHAR:   sprintf (arg, " <c>"); break;
    case TSTRING: sprintf (arg, " <s>"); break;
    case TDOUBLE: sprintf (arg, " <x>"); break;
    case TINT:    sprintf (arg, " <n>"); break;
    }

    sprintf (helpline,"  -%c%s, --%s%s:", opt [i].shortopt, arg,
	    opt [i].longopt,  arg);

    printf ("%-40s %s\n", helpline, opt [i].help);
  }
}
