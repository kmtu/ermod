#!/bin/bash

SCRIPT_WHERE=$0
SCRIPT_DIR=${SCRIPT_WHERE%/*}

failwith() {
    echo $@ 1>&2
    exit 1
}

args=$#
if (( args <= 3 )); then
    failwith "Usage: ./postprocess_namd.sh (psf) (dcd) (xst) (namd log)"
fi

PSF=$1
DCD=$2
XST=$3
LOG=$4

DCDDIR=${DCD%/*}

# check parameter consistency


# check file existence
if [[ ! -e $PSF ]]; then
    failwith "Error: PSF file does not exist!"
fi

if [[ ! -e $DCD ]]; then
    failwith "Error: DCD file does not exist!"
fi

if [[ ! -e $XST ]]; then
    failwith "Error: XST file does not exist!"
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

grep "^Info:" $LOG | perl <(cat - <<'EOF'

$temp = 300; # as default
$constant = 1;
$switchlj = 0;
$is_periodic = 0;
$coulombtype = 0;
while(<>){
# Info: BERENDSEN PRESSURE COUPLING ACTIVE
# Info: LANGEVIN PISTON PRESSURE CONTROL ACTIVE
  if(m/PRESSURE (COUPLING|CONTROL) ACTIVE/){ $constant = 2; }

# Info: LANGEVIN DYNAMICS ACTIVE
# Info: LANGEVIN TEMPERATURE   300
  if(m/LANGEVIN TEMPERATURE\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $temp = $1; }
  if(m/BERENDSEN TEMPERATURE\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $temp = $1; }

# Info: SWITCHING ACTIVE
# Info: SWITCHING ON           9
# Info: SWITCHING OFF          10
  if(m/SWITCHING OFF\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $ljcut = $1; }
  if(m/SWITCHING ON\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/){ $switchlj = $1; }

  $elecut = $ljcut;

# Info: PERIODIC CELL BASIS 1  64 0 0
  if(m/PERIODIC CELL BASIS/){ $is_periodic = 1; $coulombtype = 1; }

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

print "&ene_param
      slttype = 1,
      estype = $constant,
      boxshp = $is_periodic,
      inptemp = $temp,
      elecut = 12.0,
      ljformat = 1,
      upljcut = $ljcut,
      lwljcut = $switchlj,
      cltype = $coulombtype,
      screen = $alpha,
      splodr = $pme_order,
      ms1max = $box1,
      ms2max = $box2,
      ms3max = $box3,
      engdiv = 1,
      block_threshold = 5.0
/
&hist
      eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
      ecdmin=-80.0e0, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
      ecdmax=1.0e11, pecore=200
/
"
EOF
) > parameters_er || failwith "Failed to generate parameter from log"

# TODO: use below to check LJ mean type
# Info: USING ARITHMETIC MEAN TO COMBINE L-J SIGMA PARAMETERS"

if [[ -e engsln.01 ]] || [[ -e engsln.tt ]]; then
    echo "Warning: previous output remains, are you running the program twice?" 1>&2
fi
ln -sf $DCD ./HISTORY
ln -sf $XST ./HISTCELL


