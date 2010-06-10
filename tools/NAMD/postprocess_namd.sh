#!/bin/bash

SCRIPT_WHERE=$0
SCRIPT_DIR=${SCRIPT_WHERE%/*}

failwith() {
    echo $@ 1>&2
    exit 1
}

pos=0
while true; do
    arg=$1
    shift
    if [[ $arg = "" ]]; then
	break
    fi
    case $arg in
	(--div=*)
	DIV=${arg##--div=}
	;;
	(--*)
	failwith "unknown argument $arg"
	;;
	(*)
	case $pos in
	    (0)
	    PSF=$arg
	    ;;
	    (1)
	    DCD=$arg
	    ;;
	    (2)
	    LOG=$arg
	    ;;
	    (3)
	    DCDVAC=$arg
	    ;;
	    (*)
	    failwith "too many arguments"
	    ;;
	esac
	(( pos = pos + 1 ))
	;;
    esac
done

if [[ $LOG = "" ]]; then
    echo "Usage: ./postprocess_namd.sh (psf) (dcd) (namd log) [--div=n]" 1>&2
    echo "       ./postprocess_namd.sh (psf) (dcd) (namd log) (dcd in vacuo) [--div=n]" 1>&2
    exit 1
fi

XST=${DCD%.dcd}.xst
DCDDIR=${DCD%/*}

# echo PSF $PSF DCD $DCD LOG $LOG DCDVAC $DCDVAC

# check parameter consistency


# check file existence
if [[ ! -e $PSF ]]; then
    failwith "Error: PSF file does not exist!"
fi

if [[ ! -e $DCD ]]; then
    failwith "Error: DCD file does not exist!"
fi

if [[ ! -e $XST ]]; then
    echo "Warning: XST file does not exist, are you running non-periodic simulation?"
fi

if [[ ! -e $LOG ]]; then 
    failwith "Error: NAMD log file does not exist!"
fi

if [[ $DCDVAC != "" ]] && [[ ! -e $DCDVAC ]]; then
    failwith "Error: DCD file in vacuo is specified but does not exist"
fi

if [[ ! -e MDinfo ]]; then
    failwith "MDinfo file does not exist; make file according to the document"
fi

if [[ ! -e SltInfo ]]; then
    failwith "SltInfo file does not exist; try generating it by gen_sltinfo"
fi

# checking condition

if [[ $DCDVAC = "" ]]; then
    # in solution, no insertion
    INS=0
else
    # in refrence solution system
    INS=1
fi

# 10-5 as default
if [[ $DIV = "" ]]; then
    if [[ $INS = 0 ]]; then
	DIV=10
    else
	DIV=5
    fi
fi

grep "^Info:" $LOG | DIV=$DIV INS=$INS perl <(cat - <<'EOF'

$temp = 300; # as default
$constant = 1;
$switchlj = 0;
$is_periodic = 0;
$coulombtype = 0;
if($ENV{INS} == 1){ $instype = 3; } else { $instype = 1; }
$numdiv=$ENV{DIV};
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
      slttype = $instype,
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
      engdiv = $numdiv,
      block_threshold = 5.0
/
&hist
      eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
      ecdmin=-20.0e0, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
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

if [[ -e $XST ]]; then
    ln -sf $XST ./HISTCELL
fi

if [[ $INS = 1 ]]; then
    ln -sf $DCDVAC ./SltConf
fi


