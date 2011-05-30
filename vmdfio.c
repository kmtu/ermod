/*
  Generic trajectory access via VMD plugin.
 */

/* FIXME: add F77_FUNC(small,CAPITAL) issues */
#include "config.h"

#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <math.h>

#include <dlfcn.h>

#include "vmdplugins/molfile_plugin.h"

/*
  Declaration of plugin access points.
*/

#ifndef INSTALL_PLUGIN_PATH
#warning INSTALL_PLUGIN_PATH is not defined
#define INSTALL_PLUGIN_PATH ""
#endif

#define ADHOC_RESOLVER "adhocresolve.list"

/* lazy */
#define MAXTYPES 1000
static int typecounts = 0;
static molfile_plugin_t* vmdpluginentries[MAXTYPES];

/* lazy too */
#define MAXPLUGINS 100
static int plugincounts = 0;
static void* ldhandles[MAXPLUGINS];

/*
  Callback function. Must be compatible with vmdplugin_register_cb. 
 */
static int register_callback(void* __vp, vmdplugin_t * entity)
{
  if(entity == NULL || entity -> type == NULL) return 0;
  vmdpluginentries[typecounts++] = (molfile_plugin_t*)entity;
  return 0;
}

static int check_directory_available(char *path)
{
  struct stat st;
  int r = stat(path, &st);

  if(r == 0 && S_ISDIR(st.st_mode))
    return 1;
  return 0;
}

static int is_executable(char* path)
{
  struct stat st;
  int r;
  r = stat(path, &st);
  if(r != 0) return 0;
  if(!(S_ISREG(st.st_mode) || S_ISLNK(st.st_mode))) {
    return 0;
  }
  
  if(access(path, X_OK) == 0) return 1;
  return 0;
}

static char* strdupcat(char* str1, char* str2)
{
  char* buf = malloc(strlen(str1) + strlen(str2) + 1);
  strcpy(buf, str1);
  strcat(buf, str2);
  return buf;
}

static char* find_plugin_paths()
{
  char *pluginpath;
  /* contract: return malloc'ed string, or return NULL if error */

  /* is ERMOD_PLUGINS present and directory? */
  if((pluginpath = getenv("ERMOD_PLUGINS")) && 
     check_directory_available(pluginpath)) {
    return strdup(pluginpath);
  }

  /* is INSTALL_PLUGIN_PATH available? */
  if(strlen(INSTALL_PLUGIN_PATH) > 0 &&
     check_directory_available(INSTALL_PLUGIN_PATH)) {
    return strdup(INSTALL_PLUGIN_PATH);
  }

  return NULL;
}


void vmdfio_init_traj_(void)
{
  char* plugindir;
  char* buf;
  
  if((plugindir = find_plugin_paths()) == NULL){
    fprintf(stderr, "Error: cannot find plugin paths, please set ERMOD_PLUGINS environment variable\n");
    exit(EXIT_FAILURE);
  }

  /* traverse directory */
  do{
    DIR* dirp;
    struct dirent *ent;
    dirp = opendir(plugindir);
    
    while((ent = readdir(dirp)) != NULL){
      /* concat paths */
      char* d1 = strdupcat(plugindir, "/");
      char* pluginpath = strdupcat(d1, ent -> d_name);
      free(d1);
      
      if(is_executable(pluginpath) && pluginpath[strlen(pluginpath)-1] != '.'){
	void* handle;
	handle = dlopen(pluginpath, RTLD_NOW | RTLD_GLOBAL);
	if(!handle){
	  fprintf(stderr, "Warning: failed to load plugin \"%s\", reason: \"%s\"\n", pluginpath, dlerror());
	  continue;
	}
	void* initptr;
	void* entptr;
	initptr = dlsym(handle, "vmdplugin_init");
	entptr = dlsym(handle, "vmdplugin_register");
	if(initptr != NULL && entptr != NULL){
	  int r;
	  r = ((int (*)(void))initptr)();
	  if(r != 0)
	    fprintf(stderr, "Warning: error while initializing %s\n", pluginpath);
	  ldhandles[plugincounts++] = handle;
	  r = ((int (*)(void*, vmdplugin_register_cb))entptr)(NULL, register_callback);
	  if(r != 0)
	    fprintf(stderr, "Warning: error while registering %s\n", pluginpath);
	}else{
	  fprintf(stderr, "Warning: failed to load entry point (plugin: %s)\n", pluginpath);
	}
      }
    }
  }while(0);
}

void vmdfio_fini_traj_(void)
{
  int i;
  for(i = 0; i < plugincounts; ++i){
    void* entptr = dlsym(ldhandles[i], "vmdplugin_fini");
    if(entptr){
      ((int (*)(void))entptr)();
    }
    dlclose(ldhandles[i]);
  }
}

typedef struct vmdpluginio_t {
  molfile_plugin_t * plugin;
  void* filehandle;
  int natoms;
} vmdpluginio;

void vmdfio_open_traj_(void **handle, char *fname, int *fnamelen, int *status)
{
  char* buf;
  size_t buflen = 8192;
  ssize_t r;
  int fp;
  char* ext;
  int i;

  buf = malloc(sizeof(char) * buflen + 1);
  strncpy(buf, fname, *fnamelen);
  buf[*fnamelen] = '\0';

  while(1){
    r = readlink(buf, buf, buflen);
    if(r == -1) break;
    buf[r] = '\0';
  }

  if(*fnamelen == strlen(buf) && strncmp(fname, buf, *fnamelen) == 0){
    /* not a symbolic link? */
    fprintf(stderr, "Error: vmdfio.c: failed to open with vmdfio_open_traj_. (filename = \"%s\".) Perhaps it's not a symbolic link?\n", buf);
    *status = -1;
    goto cleanup;
  }
  
  /* 
     select plugin to use.
     use filename_extension to select
  */
  *status = -1;
  *handle = NULL;
  {
    /* get filename extension */
    char* lastdot = strrchr(buf, '.');

    if(lastdot == NULL) goto cleanup;
    ext = lastdot + 1;
  }
  fprintf(stderr, "Opening: %s...", buf);

  for(i = 0; i < typecounts; ++i){
    int extlen = strlen(ext);
    molfile_plugin_t *p = vmdpluginentries[i];
    char* plugin_supportext;
    char* tokptr;
    char* ptr;

    if(strcmp(p -> type, MOLFILE_PLUGIN_TYPE) != 0) continue; 
    plugin_supportext = strdup(p -> filename_extension);

    while(1){
      ptr = strtok_r(plugin_supportext, ",", &tokptr);
      if(ptr == NULL) break;
      plugin_supportext = NULL; /* in the next subsequent call NULL is passed */

      if(strncmp(ptr, ext, extlen) == 0 && 
	 (ptr[extlen] == '\0' || ptr[extlen] == ',')){
	void* fh;
	vmdpluginio* pp = malloc(sizeof(vmdpluginio));
	pp -> plugin = p;

	/* Found, open with this plugin */
	pp -> natoms = MOLFILE_NUMATOMS_UNKNOWN;
	fh = (p -> open_file_read)(buf, ext, &(pp->natoms));
	if(fh == NULL){
	  fprintf(stderr, "Error while opening %s\n", buf);
	  exit(1);
	}
	pp -> filehandle = fh;
	fprintf(stderr, "OK, successfully opened with \"%s\"\n", p -> prettyname);

	*handle = pp;
	break;
      }
    }
    free(plugin_supportext);
    if(*handle){
      *status = 0;
      break;
    }
  }

  if(*status == -1){
    fprintf(stderr, "Error: could not find appropriate plugins\n");
  }


 cleanup:
  free(buf);
  return;
}

void vmdfio_read_traj_step_(void **handle, double* xout, double* box, int *natoms_aux, int *status)
{
  vmdpluginio *p = *handle;
  molfile_plugin_t *plugin = p -> plugin;
  molfile_timestep_t snapshot;
  int natoms = p -> natoms;
  float* buf;
  int r;
  static int firstcall = 1;

  if(natoms == MOLFILE_NUMATOMS_UNKNOWN){
    /* not determined from the trajectory */
    natoms = *natoms_aux;
  }else{
    /* check integrity */
    if(natoms != p -> natoms && firstcall){
      fprintf(stderr, "Warning: # of atoms in trajectory does not match with # of atoms in configurations\n");
      firstcall = 0;
    }
  }

  buf = malloc(sizeof(float) * 3 * natoms);
  snapshot.coords = buf;

  r = (plugin -> read_next_timestep)(p -> filehandle, natoms, &snapshot);

  do{
    int i;
    double x, y, u, v, w;

    if(r == MOLFILE_EOF){
      *status = -1;
      break;
    }

    /* both are the same format */
    for(i = 0; i < natoms * 3; ++i)
      xout[i] = buf[i];

    /*
      ~a = (1, 0, 0)
      ~b = (x, y, 0)
      ~c = (u, v, w)
      ~a.~b = x = cos gamma
      |~a*~b| = y = sin gamma
      ~a.~c = u = cos beta
      ~b.~c = xu + yv = cos alpha
    */

    x = cos(snapshot.gamma * M_PI / 180.0);
    y = sin(snapshot.gamma * M_PI / 180.0);
    u = cos(snapshot.beta * M_PI / 180.0);
    v = (cos(snapshot.alpha * M_PI / 180.0) - x * u) / y; /* FIXME: potential underflow risk */
    w = sqrt(1 - u * u - v * v);  /* FIXME: same above */

    box[0] = snapshot.A;
    box[3] = snapshot.B * x; 
    box[4] = snapshot.B * y;
    box[6] = snapshot.C * u; 
    box[7] = snapshot.C * v;
    box[8] = snapshot.C * w;
    
    box[1] = 0.; box[2] = 0.; box[5] = 0;

    *status = 0;
  }while(0);

  free(buf);
  return;
}

void vmdfio_close_traj_(void **handle)
{
  vmdpluginio *p = *handle;
  molfile_plugin_t *plugin = p -> plugin;
  (plugin -> close_file_read)(p -> filehandle);
  free(*handle);
}
