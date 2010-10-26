

/* FIXME: add F77_FUNC(small,CAPITAL) issues */

#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

// Gromacs-y files. NEED -I to include/gromacs directory, not gromacs directory.
#include <statutil.h>
#include <smalloc.h>

static const char *filenames[] = { "HISTORY", "SltConf" };

typedef struct gmxfileio_t {
  int fp;
  int natoms;
  rvec *x;
} gmxfileio;

/*
  This may not work on some environment;
  There are no portable ways to pass characters b/w C and Fortran...
 */
void open_gmtraj_(void **handle, int *which, int *status)
{
  char* buf;
  size_t buflen = 8192;
  ssize_t r;
  int fp;
  gmxfileio *fh;

  buf = malloc(sizeof(char) * buflen + 1);
  r = readlink(filenames[*which], buf, buflen);

  if(r == -1){
    /* Error cases */
    if(errno == EINVAL) {
      /* not a symbolic link? */
      fprintf(stderr, "Warning: gromacsfio.c failed to open with open_traj. (Perhaps it's not a symbolic link?)\nFalling back to standard g96 I/O routine.\n");
    }
    *status = -1;

    goto cleanup;
  }
  
  buf[r] = '\0';

  { 
    real t;
    rvec *x;
    matrix box;
    int natoms;

    natoms = read_first_x(&fp, buf, &t, &x, box);
    if(natoms == 0) {
      *status = -1;
      goto cleanup;
    }
    rewind_trj(fp);

    snew(fh, 1);
    fh -> fp = fp;
    fh -> natoms = natoms;
    fh -> x = x;
  }

  *handle = fh;
  *status = 0;

  // close_trn(fp);

 cleanup:
  free(buf);
  return;
}

void read_gmtraj_step_(void **handle, double* x, double* box, int *status)
{
  gmxfileio *fh = (gmxfileio*)*handle;

  {
    real t;
    matrix boxtmp;
    bool cont;
    int i;

    cont = read_next_x(fh -> fp, &t, fh -> natoms,
		       fh -> x, boxtmp);
    if(cont == FALSE){
      fprintf(stderr, "No more frames in trajectory (Check MDinfo file)\n");
      *status = -1;
      return;
    }

    for(i = 0; i < fh -> natoms; ++i) {
      x[i * 3 + 0] = fh -> x[i][0];
      x[i * 3 + 1] = fh -> x[i][1];
      x[i * 3 + 2] = fh -> x[i][2];
    }

    for(i = 0; i < 3; ++i) {
      box[i * 3 + 0] = boxtmp[i][0];
      box[i * 3 + 1] = boxtmp[i][1];
      box[i * 3 + 2] = boxtmp[i][2];
    }

    *status = 0;
  }
  
  return;
}

void close_gmtraj_(void **handle)
{
  gmxfileio *fh = *handle;
  close_trj(fh -> fp);
  sfree(*handle);
}
