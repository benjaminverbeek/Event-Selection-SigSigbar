#!/bin/bash

file_exist(){
		if [ ! -f $1 ] ; then
				echo "The file : $1 is not exist, create..."
				return 1
		fi
}


current=`pwd`
# path="/ustcfs/BES3User/2021/benjamin/runarea/rec/2009/rec/"
path="/ustcfs/bes3data/708-1/jpsi/round02/mc/"
dstpath="outfile"
jobpath=$current/cmtfile
jobfile="jobOptions_ana_sigmasigmabar.txt"
prefix="SigmaSigmabar"
dstfilepath=$current/$dstpath/dst
tmpfilepath=$current/$dstpath/tmp
rootfilepath=$current/$dstpath/root
# njob=99
njob=470 # there are 4700 .dst-files in the folder

rdmseed=5009

file_exist data.txt
if test $? -ne 0; then
	find $path -name "*.dst" | sort > data.txt
fi
mkdir -p $dstfilepath $tmpfilepath $rootfilepath $jobpath tmp
nline=`awk 'END{print NR}' data.txt`
lcount=`expr $nline / $njob`
split -l $lcount data.txt -d -a 4 tmp/$prefix

njobcout=0
for file in tmp/*
do
		echo $file
		opt_first=true
		input=''
		for line in `cat $file`
		do
			if $opt_first; then
				opt_first=false
				input="\n\"$line\""
			else 
				input="$input, \n\"$line\""
			fi
		done
		njobcout=`expr $njobcout + 1`
		outfile=${prefix}_$(printf %04d ${njobcout}).txt
		outdst=$dstfilepath/${prefix}_$(printf %04d ${njobcout}).dst
		outtmp=$tmpfilepath/${prefix}_$(printf %04d ${njobcout}).tmp
		outroot=$rootfilepath/${prefix}_$(printf %04d ${njobcout})
		cat $jobfile | sed "s#INPUT#$input#g" | sed "s#DSTFILE#$outdst#g" | sed "s#TMPFILE#$outtmp#g" | sed "s#OUTPUT#$outroot#g" | sed "s#RANDOM#$rdmseed#g" > $jobpath/$outfile
		rdmseed=`expr $rdmseed + 1`
done

rm -rf tmp
rm -rf data.txt



if false
then
nline=`awk 'END{print NR}' data.txt`
ntmp=`expr $nline / $njob`
ntmp=`expr $ntmp + 1`
ntmp=1

ncout=0
njobcout=0
opt_first=true
input=''
mkdir -p $dstfilepath $tmpfilepath $rootfilepath $jobpath
for line in `cat data.txt`
do	
		if test $[ncout]  -eq $[ntmp]
		then
				njobcout=`expr $njobcout + 1`
				outfile=${prefix}_$(printf %04d ${njobcout}).txt
				outdst=$dstfilepath/${prefix}_$(printf %04d ${njobcout}).dst
				outtmp=$tmpfilepath/${prefix}_$(printf %04d ${njobcout}).tmp
				outroot=$rootfilepath/${prefix}_$(printf %04d ${njobcout})
				input="$input, \n\"$line\""
				cat $jobfile | sed "s#INPUT#$input#g" | sed "s#DSTFILE#$outdst#g" | sed "s#TMPFILE#$outtmp#g" | sed "s#OUTPUT#$outroot#g" > $jobpath/$outfile
				echo $jobpath/$outfile
			#	cat $jobfile  > $dstpath/$outfile
				ncout=0
				input=''
				opt_first=true
		else
				ncout=`expr $ncout + 1`
				if $opt_first; then
						opt_first=false
						input="\n\"$line\""
				else
						input="$input, \n\"$line\""
				fi
		fi
done
fi





