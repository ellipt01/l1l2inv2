#!/usr/bin/env bash

#### function definitions ####

# display usage and exit
usage_exit() {
	echo ""
	echo "PROGAM: ${0##*/}"
	echo ""
	echo "DESGRIPTION:"
	echo "       This script derives interpolated L-curve and its curvature"
	echo "       with the aid of smooth b-spline"
	echo ""
	echo "USAGE: ${0##*/}"
	echo "[optional]"
	echo "       -f <regression info filename: default is regression_info.data>"
	echo "       -a <alpha: default is 0.9>"
	echo "       -b <regularization parameter for smooth spline; default is 0.01>"
	echo "          if -q option is not specified"
	echo "          (i.e., second-step inversion is not performed),"
	echo "          this option is ignored"
	echo "       -c (return only values:"
	echo "          optimal lambda, RSS, solution norm, max_curvature,"
	echo "          and iteration number of optimal lambda)"
	echo "       -v (verbose mode)"
	echo "       -h (show this message and exit)"
	exit 1
}

dir=`dirname ${BASH_SOURCE}`

###
alpha=0.9
beta=0.01
fn="regression_info1.data"
while getopts "f:a:b:cvh" OPT; do
	case $OPT in
		f)  fn=$OPTARG ;;
		a)  alpha=$OPTARG ;;
		b)  beta=$OPTARG ;;
		c)  only_values=1 ;;
		v)  VERBOSE=1 ;;
		h)  usage_exit ;;
		/?) usage_exit ;;
	esac
done

if [ ! -z $VERBOSE ]; then
	echo "$dir"/lcurve_interp -f $fn -a $alpha -b $beta
fi
"$dir"/lcurve_interp -f $fn -a $alpha -b $beta


if [ -e splined.data ]; then
	mv -f splined.data splined1.data
fi

rm -f curvature.data
mv cv.data curvature.data

# find optimal lambda
res=`cat curvature.data | gawk 'BEGIN{max=0;l=0;j=0}\
	NR>1{if($2>max){max=$2;l=$1;j=NR+2}}\
	END{print l,max,j}'`
lambda=`echo $res | gawk '{print $1}'`
curv=`echo $res | gawk '{print $2}'`
lj=`echo $res | gawk '{print $3}'`

if [ `echo "$alpha > 0." | bc` == 1 ]; then
	# L1L2 norm regularization
	lval=`gawk 'BEGIN{print '$lambda'*'$alpha'}'`
	ltype=1
elif [ `echo "$alpha == 0." | bc` == 1 ]; then
	# L2 norm regularization
	lval=$lambda
	ltype=2
fi

lx=`echo $lval | gawk '{print round(log($1)/log(10))}\
	function round(x){return (x>0)?int(x+0.5):int(x-0.5)}'`
eps=`gawk 'BEGIN{print 10^('$lx'-2)}'`

if [ $ltype -eq 1 ]; then
	res=`cat $fn\
		| gawk 'NR>1{if(abs($4-'$lval')<'$eps'){print $3,$1,NR-2;exit}}\
		function abs(x){return (x<0)?-x:x}'`
elif [ $ltype -eq 2 ]; then
	res=`cat $fn\
		| gawk 'NR>1{if(abs($5-'$lval')<'$eps'){print $3,$2,NR-2;exit}}\
		function abs(x){return (x<0)?-x:x}'`
fi

rss=`echo $res | gawk '{print $1}'`
nrm=`echo $res | gawk '{print $2}'`
ljopt=`echo $res | gawk '{print $3}'`

if [ -z $only_values ]; then
	echo "lambda = " $lambda
	echo "rss = " $rss
	echo "nrm = " $nrm
	echo "curvature = " $curv
	echo "optimal: lambda["$ljopt"] =" $lambda
else
	echo $lambda $rss $nrm $curv $ljopt
fi

