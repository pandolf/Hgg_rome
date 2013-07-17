#!/bin/csh
# $Id: preparelist_eth.csh,v 1.3 2011/07/29 09:16:48 meridian Exp $

if( $#argv<3  ) then
  echo "usage:  prepareList.csh  <inputfile> <listname> <location>  [run if 1]"
  exit 0
endif

set run = 0
if( $#argv>2 ) then
  set run = $4
endif

set infile = $1

set listname = $2

set location = $3

# num of files per list file
set filexlist  = 20

set prepend=""

if ( $location == "cern" ) then
  set prepend=""
else if ( $location == "xrootd" ) then
  set prepend="root://pccmsrm23.cern.ch:1094/"
else if ( $location == "eth" ) then
  set prepend="dcap://t3se01.psi.ch/"
endif 

#set prepend="dcap://cmsrm-se01.roma1.infn.it"

 
#set files = ( `/bin/ls -1 $indir` )
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Working on dataset $listname" 
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

if ($location != "xrootd" ) then
    set files = ( `cat $infile | grep $listname` )
else
    set files = ( `cat $infile | grep $listname | sed -e "s%/xrootdfs%%g"` )
endif


set tmpfile = /tmp/tmpfilelist.${listname}
rm -Rf $tmpfile
touch $tmpfile

foreach i ( $files )
    echo ${prepend}$i >> ${tmpfile}
end

set unifile =  /tmp/uniq.${listname}.files
# Removing duplicates and using the last file (according to order in the inputfile) for a given jobId
awk 'BEGIN{FS="output_"}{ split($2,jobId,"_"); sum[jobId[1]]++;line[jobId[1]]=$0}END{for(i in sum) print line[i];}' ${tmpfile} | sort -n >! ${unifile}.tmp

#rm -Rf  ${unifile}.tmp
#sort $tmpfile | uniq -w 4 | awk -F: '{print $2":"$3}' > ${unifile}.tmp

set suffixlen = 2

split -l $filexlist -d -a $suffixlen  ${unifile}.tmp  ${listname}_

set uniqfiles=`wc -l ${unifile}.tmp`

echo "# of root files in directory: $#files"
echo "# of uniq files in directory: ${uniqfiles}"
echo "# of files per list: $filexlist"

foreach i ( ${listname}_?? )
  mv $i ${i}.list
  echo "new list:   ${i}.list"
end



if( $run != 1 ) then
  rm -rf ${listname}_*.list
endif


#rm -Rf $tmpfile

