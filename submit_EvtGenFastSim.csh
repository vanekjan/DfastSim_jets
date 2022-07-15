#!/bin/csh

rm -r LocalLibraries.package
rm LocalLibraries.zip

set prodId=`date +%F_%H-%M`

mkdir ./myOutput/${prodId}
mkdir ./jobs/log/${prodId}
mkdir ./jobs/err/${prodId}
mkdir ./jobs/submit/${prodId}
mkdir ./jobs/report/${prodId}
mkdir ./jobs/csh/${prodId}
mkdir ./jobs/list/${prodId}

star-submit-template -template submitToyMcEvtGen.xml -entities productionId=${prodId}

