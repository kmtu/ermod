#!/bin/bash

# rewrite the following variables according to your environment
PROGDIR=.. # where program is (can be either absolute or relative path)
PROGCONF=$PROGDIR

SCRIPT_WHERE=$0
SCRIPT_DIR=${SCRIPT_WHERE%/*}

failwith() {
    echo $@ 1>&2
    exit 1
}

args=$#
if (( args <= 2 )); then
    failwith "Usage: ./postprocess_namd.sh (psf) (dcd) (namd log)"
fi

PSF=$1
DCD=$2
LOG=$3

DCDDIR=${DCD%/*}

# check environment consistency
if [[ -e $PROGDIR/enganal.f ]]; then
    failwith "Error: PROGDIR is not program directory; edit PROGDIR variable in $0"
fi

if [[ -e $PROGDIR/Makefile ]]; then
    failwith "Error: energy representation program is not configured; try ./configure on $PWD/$PROGDIR"
fi

# check parameter consistency


# check file existence
if [[ ! -e $PSF ]]; then
    failwith "Error: PSF file does not exist!"
fi

if [[ ! -e $DCD ]]; then
    failwith "Error: DCD file does not exist!"
fi

if [[ ! -e $LOG ]]; then 
    failwith "Error: NAMD log file does not exist!"
fi

if [[ ! -e MDinfo ]]; then
    failwith "MDinfo file does not exist; make file according to the document"
fi

if [[ ! -e SltInfo ]]; then
    failwith "SltInfo file does not exist; try generating it by inpfile.f"
fi

grep "^Info:" $LOG | perl -e <(cat - <<'EOF'

while(<>){
# Info: BERENDSEN PRESSURE COUPLING ACTIVE
# Info: LANGEVIN PISTON PRESSURE CONTROL ACTIVE
  if(m/PRESSURE (COUPLING|CONTROL) ACTIVE/){ $estype = 2; }else{ $estype = 1; }
  $temp = 300; # as default
# Info: LANGEVIN DYNAMICS ACTIVE
# Info: LANGEVIN TEMPERATURE   300
  if(m/LANGEVIN TEMPERATURE\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $temp = $1; }
  if(m/BERENDSEN TEMPERATURE\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $temp = $1; }

# Info: SWITCHING ACTIVE
# Info: SWITCHING ON           9
# Info: SWITCHING OFF          10
  if(m/SWITCHING OFF\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $ljcut = $1; }
  if(m/SWITCHING ON\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $switchlj = $1; }else{ $switchlj = $ljcut }

  $elecut = $ljcut

# Info: PERIODIC CELL BASIS 1  64 0 0
  if(m/PERIODIC CELL BASIS/){ $coulombtype = 1; } else { $coulombtype = 0; }

# Info: PARTICLE MESH EWALD (PME) ACTIVE
# Info: PME TOLERANCE               1e-06"
# Info: PME EWALD COEFFICIENT       0.3"
# Info: PME INTERPOLATION ORDER     4"
# Info: PME GRID DIMENSIONS         64 64 128"
  if(m/\(PME\) ACTIVE/){ $coulombtype = 2; }
  if(m/PME EWALD COEFFICIENT\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $alpha = $1; }
  if(m/PME INTERPOLATION ORDER\s+(\d+)/){ $pme_order = $1; }
  if(m/PME GRID DIMENSIONS\s+(\d+) (\d+) (\d+)/){ $box1 = $1; $box2 = $2; $box3 = $3; }
# FIXME: how to obtain parameters for Ewald version?

}

print "c
      slttype = 1
      estype = $constant
      boxshp = $is_periodic
      inptemp = $temp
      elecut = 12.0
      upljcut = $ljcut
      lwljcut = $switchlj
      cltype = $coulombtype
      screen = $alpha
      splodr = $pme_order
      ms1max = $box1
      ms2max = $box2
      ms3max = $box3
      engdiv = 1
      eclbin=5.0e-2 ; ecfbin=2.0e-3 ; ec0bin=2.0e-4 ; finfac=10.0e0
      ecdmin=-80.0e0 ; ecfmns=-0.20e0 ; ecdcen=0.0e0 ; eccore=20.0e0
      ecdmax=1.0e11 ; pecore=200
      box_threashold = 5.0
"
EOF
) > $PROGDIR/param_eng

# TODO: use below to check LJ mean type
# Info: USING ARITHMETIC MEAN TO COMBINE L-J SIGMA PARAMETERS"

pushd $PROGDIR
make || failwiwth "Error compling program"
popd

if [[ ( -e engsln.01 ) -o ( -e engsln.tt ) ]]; then
    echo "Warning: previous output remains, are you running the program twice?" 1>&2
fi
ln -sf $DCD ./HISTORY
$PROGDIR/er_namd || failwith "Error executing main"


