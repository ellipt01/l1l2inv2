#!/usr/bin/env bash


#### function definitions ####

# display usage and exit
usage_exit() {
	echo ""
	echo "PROGAM: ${0##*/}"
	echo ""
	echo "DESGRIPTION:"
	echo "       This script performs L1-L2 norm magnetic inversion"
	echo "       based on CDA with warm-start by calling \"l1l2inv\","
	echo "       as well as estimates the most optimal regularization parameter"
	echo "       by estimating the curvature of the L-curve."
	echo ""
	echo "USAGE: ${0##*/}"
	echo "       -a <alpha>"
	echo "[optional]"
	echo "       -w <lambda_min:dlambda; default is -1:0.1>"
	echo "          inversion is performed for a sequential lambda"
	echo "          starting from log10(lambda_max) down to log10(lambda_min)"
	echo "          with an interval of log10(dlambda),"
	echo "          where, lambda_max is automatically estimated"
	echo "          by main inversion program \"l1l2inv\"" 
	echo "       -t <tolerance; default=1.e-5>"
	echo "       -n <bounds of solution lower:upper; default is no bounds>"
	echo "       -s <parameter settings file; default is ./settings>"
	echo "       -b <regularization parameter for smooth spline; default is 0.01>"
	echo "          if -q option is not specified"
	echo "          (i.e., second-step inversion is not performed),"
	echo "          this option is ignored"
	echo "       -g <0(false) or 1(true): stretch the grid cells"
	echo "           at the edge of the model space outward, default is 1>"
	echo "       -p (perform CDA using parallel computing; default is none)"
	echo "       -q (perform second-step inversion using most opt-lambda;"
	echo "           default is none)"
	echo "       -c (use stochastic CDA instead of cyclic CDA; default is none)"
	echo "       -v (verbose mode)"
	echo "       -h (show this message and exit)"
	exit 1
}

# set inline options for core inversion program "l1l2inv"
setopts() {
	OPTS=""

	if [ ! -z "$WEIGHTS" ]; then
		OPTS="-d $WEIGHTS $OPTS"
	fi

	if [ ! -z "$BOUNDS" ]; then
		OPTS="-n $BOUNDS $OPTS"
	fi

	if [ ! -z "$GRID" ]; then
		OPTS="-g $GRID $OPTS"
	fi

	if [ ! -z $STOCHASTIC ]; then
		OPTS="$OPTS -c"
	fi

	if [ ! -z $PARALLEL ]; then
		OPTS="$OPTS -p"
	fi

	if [ ! -z $OUTPUT_VECTORS ]; then
		OPTS="$OPTS -o"
	fi

	if [ ! -z $OUTPUT_WEIGHTED ]; then
		OPTS="$OPTS -u"
	fi

	if [ ! -z $ACTUALMAGPENALTY ]; then
		OPTS="$OPTS -k"
	fi

	if [ ! -z "$VERBOSE" ]; then
		OPTS="$OPTS -v"
	fi
	echo "$OPTS"
}


#### starting script ####

if [ -z "$1" ]; then
	usage_exit
fi

dir=`dirname ${BASH_SOURCE}`

echo $0 $@
echo ""

RANGE="-1:0.1"
TOL=1.e-5
SFILE="./settings"
BETA=0.01
TYPE=1 # L1L2

while getopts "r:d:a:w:t:n:s:b:g:kpqcouvh" OPT; do
	case "$OPT" in
		r)  TYPE=$OPTARG ;;
		d)  WEIGHTS=$OPTARG ;;
		k)  ACTUALMAGPENALTY=1 ;;
		a)  ALPHA=$OPTARG ;;
		w)  RANGE=$OPTARG ;;
		t)  TOL=$OPTARG ;;
		n)  BOUNDS=$OPTARG ;;
		s)  SFILE=$OPTARG ;;
		b)  BETA=$OPTARG ;;
		g)  GRID=$OPTARG ;;
		p)  PARALLEL=1 ;;
		q)  SPLINE=1 ;;
		c)  STOCHASTIC=1 ;;
		o)  OUTPUT_VECTORS=1 ;;
		u)  OUTPUT_WEIGHTED=1 ;;
		v)  VERBOSE=1 ;;
		h)  usage_exit ;;
		/?) usage_exit ;;
	esac
done 

if [ -z $ALPHA ]; then
	echo "ERROR: please specify alpha as -a <alpha>"
	usage_exit
fi

if [ ! -e input.data ]; then
	echo "ERROR: ${0##*/}: file \"input.data\" not exists. Abort!"
	exit 1
fi

if [ -z "$SFILE" -a ! -e "settings" ]; then
	echo "ERROR: ${0##*/}: file \"settings\" not exists."
	echo "if you use your own setting file, please specify as -s <filename>."
	exit 1
fi

if [ ! -e "$dir"/l1l2inv ]; then
	echo "ERROR: ${0##*/}: inversion program \"l1l2inv\" is not found"
	echo "in the directory where this script exists."
	echo "please compile and install programs correctly"
	echo "by performing make, make install, etc."
	exit 1
fi

OPTS=`setopts`
echo -n >| optimal_info.data

#############################################################
#    #### first step inversion ####
# optimization with respect to a sequential lambda
# starting from log10(lambda_max)
# down to a log10(lambda_min)
# with an interval of log10(dlambda)
#
# lambda_max: automatically estimated
# lambda_min, dlambda: specified by user
#############################################################

echo "$dir"/l1l2inv -a $ALPHA -w "$RANGE" -t $TOL -s $SFILE ${OPTS}
"$dir"/l1l2inv -a $ALPHA -r $TYPE -w "$RANGE" -t $TOL -s $SFILE ${OPTS}

if [ $? -ne 0 ]; then
	echo
	echo "ERROR: program \"l1l2inv\" did not finish successfully. script abort!"
	echo
	exit 1
fi

# evaluate curvature of the L-curve
"$dir"/curvature.sh -f regression_info.data -a $ALPHA -b $BETA >> optimal_info.data
cat optimal_info.data

# store prev inversion results
if [ -e "regression_info.data" ]; then
	mv regression_info.data regression_info1.data
fi
if [ -e "beta_path.data" ]; then
	mv beta_path.data beta_path1.data
fi
if [ -e "curvature.data" ]; then
	mv curvature.data curvature1.data
fi
if [ -e "lcurve_interp.data" ]; then
	mv lcurve_interp.data lcurve_interp1.data
fi


### following codes are executed if -q option is specified ###

#############################################################
#    #### second step inversion (optional) ####
# optimization for most optimal lambda
# with an aid of the smooth-spline(b-spline)
#############################################################

if [ ! -z $SPLINE ]; then

	# eval optimal lambda using b-spline
	lambda=`"$dir"/optimal_lambda -f curvature1.data -d 0.001`
	echo "most optimal lambda =" $lambda >> optimal_info.data
	echo "most optimal lambda =" $lambda

	if [ `echo "$ALPHA > 0." | bc` == 1 ]; then
		# L1L2 norm regularization
		lval=`gawk 'BEGIN{print '$lambda'*'$ALPHA'}'`
		ltype=1
	elif [ `echo "$ALPHA == 0." | bc` == 1 ]; then
		# L2 norm regularization
		lval=$lambda
		ltype=2
	fi

	if [ $ltype -eq 1 ]; then
		res=`cat regression_info1.data\
			| gawk '{if(NR>1){\
				if(NR==2){\
					l=$4;k=0\
				}else{\
					if($4<'$lval'){print k,l;exit}else{l=$4;k=k+1}\
				}\
			}}'`
	elif [ $ltype -eq 2 ]; then
		res=`cat regression_info1.data\
			| gawk '{if(NR>1){\
				if(NR==2){\
					l=$5;k=0\
				}else{\
					if($5<'$lval'){print k,l;exit}else{l=$5;k=k+1}\
				}\
			}}'`
	fi
	i=`echo $res | gawk '{print $1}'`
	"$dir"/extract -f beta_path1.data -i $i -o 1 -s $SFILE >| beta"$i".data

	lmax=`echo $res | gawk '{print $2}'`
	if [ "$ltype" -eq 1 ]; then
		lmax=`gawk 'BEGIN{print '$lmax'/'$ALPHA'}'`
	fi

	# perform second-step inversion
	lmin=$lambda
	llmax=`gawk 'BEGIN{print log('$lmax')/log(10)}'`
	llmin=`gawk 'BEGIN{print log('$lmin')/log(10)}'`
	dll=`gawk 'BEGIN{print (('$llmax')-('$llmin'))/5}'`
	RANGE=`gawk 'BEGIN{printf("%.8e:%.8e:%.8e",'$llmin','$dll','$llmax')}'`
	echo "$dir"/l1l2inv -a $ALPHA -w "$RANGE" -t $TOL -s $SFILE -b beta"$i".data ${OPTS}
	"$dir"/l1l2inv -a $ALPHA -r $TYPE -w "$RANGE" -t $TOL -s $SFILE -b beta"$i".data ${OPTS}

	if [ $? -ne 0 ]; then
		echo
		echo "ERROR: program \"l1l2inv\" did not finish successfully. script abort!"
		echo
		exit 1
	fi
fi

exit 0
