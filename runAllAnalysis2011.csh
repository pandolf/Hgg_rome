#!/bin/tcsh

#set data_json = "`pwd`/jsonFiles/Cert_160404-180252_7TeV_All2011_v3.txt"
#set data_json = "`pwd`/jsonFiles/latest_prompt_2012.json"
set data_json = "`pwd`/jsonFiles/Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2.txt"
set puweight_41x = "`pwd`/mc_41x_PUweight.root"
set puweight_42x = "`pwd`/RUN2011_0100_73500.weights.root"
#set puweight_52x = "`pwd`/Summer12-190456-194076_Prompt_RUN2012.68300.true.weights.root"
#set puweight_52x = "`pwd`/Summer12-rereco_prompt_2012.69000.observed.weights.root"
#set puweight_52x = "`pwd`/Summer12-prompt08Jun.71000.observed.weights.root"
set puweight_52x = "`pwd`/nPU-Summer12_DD3.190456-201678-29Jun_Prompt.observed.root"
set ptweightfile_template = "`pwd`/kfactors/Kfactors_MASSVALUE_AllScales.root"

set location = "eth"
set version = "v1"
set run = 0

if($#argv == 0 || $#argv < 3 || $#argv > 7 ) then
  echo "usage:  runAllAnalysis2011.csh  <location> <version> <run if 1> <jsonfile> <pureweight> <ptweight> <energy correction/smearing> "
  echo "        locations: cern roma eth"
  echo "        version: version string for redntp"
  echo "        run: default=0  set to 1 to execute"
  echo "        jsonfile: optional json to select good RUN_LS used for data (full path is required)"
  echo "        pu weight for MC: default=-1  set to 1 to store them"
  echo "        pt weight for MC (HiggsPt): default=-1  set to 1 to store them"
  exit -1
endif

set location = $1
echo "location : $location "

set version = $2
echo "version : $version "

set run = $3
echo "run : $run "

if ($#argv > 3) then
  set data_json = $4
  echo "json : ${data_json}"
endif 

set puweight = -1
if ($#argv > 4) then
  set puweight = $5
  echo "pu weight: ${puweight}"
endif 

set ptweight = -1
if ($#argv > 5) then
  set ptweight = $6
  echo "pt weight: ${ptweight}"
endif 

set energyCorrection = -1
if ($#argv > 6) then
  set energyCorrection = $7
  set energyCorrectionName = `basename ${energyCorrection} .dat`
  echo "energyCorrection: ${energyCorrection}"
endif 

#foreach class ( 53x_globe_synch )
foreach class ( 53x_globe_synch_data )
#foreach class ( 53xv2 52xv5 )
#foreach class ( 53xv2 52xv5 53xv2_data )
    foreach preseltype ( preselectionMVA )
#    foreach preseltype ( cicpfloose )
#   foreach preseltype ( preselectionCS cicpfloose preselectionMVA cicpfloosenoeleveto )
#    foreach preseltype ( cicpfloose preselectionMVA  ) 
#    foreach preseltype ( cicpfloosenoeleveto ) 
	if ( "`echo ${class} | grep data`XXX" != "XXX" ) then
	    set command="./makeRedNtp.csh list.${class}/ redntp.${class}.${preseltype}.${energyCorrectionName}.${version} ${preseltype} ${location} ${run} $data_json -1 -1 ${energyCorrection}"
	else 
	    if ( $puweight !=  -1 ) then
		if ( "`echo ${class} | grep 41x`XXX" != "XXX" ) then
		    set puweightFile = ${puweight_41x}
		else if ( "`echo ${class} | grep 42x`XXX" != "XXX" ) then
			set puweightFile = ${puweight_42x}			
		else if ( "`echo ${class} | grep 52x`XXX" != "XXX" ) then
			set puweightFile = ${puweight_52x}
		else if ( "`echo ${class} | grep 53x`XXX" != "XXX" ) then
			set puweightFile = ${puweight_52x}
		endif
	    else
		set puweightFile = -1
	    endif

	    if ( $ptweight !=  -1 ) then
		set ptweightFile = ${ptweightfile_template}
	    else
		set ptweightFile = -1
	    endif
	    set command="./makeRedNtp.csh list.${class}/ redntp.${class}.${preseltype}.${energyCorrectionName}.${version} ${preseltype} ${location} ${run} -1 ${puweightFile} ${ptweightFile} ${energyCorrection}.MC"
	endif
	echo ${command}
	if ( $run == 1 ) then
	   ${command}
	endif
    end
end
