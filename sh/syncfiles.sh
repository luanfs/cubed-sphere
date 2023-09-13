#!/bin/bash

# Creates a backup sync of main model files
# P. Peixoto - Jul 2012
# modified by Luan Santos - 2022

#output="fvcs$version.tar.bz2"
output="fvcs.tar.bz2"
bkpdir="fvcs" 

#-------------------------------------------------------------------------------------------------------
# remote host 1 - ime.usp.br
user_remote_host1="luansantos"
remote_host1="brucutu.ime.usp.br"
remote_host1_dir="/var/tmp/lfs"
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# remote host 2 - ybytu
user_remote_host2="luansantos"
remote_host2="ybytu.ime.usp.br"
remote_host2_dir="/home/luansantos/doc/fvcs"
#-------------------------------------------------------------------------------------------------------

 
#-------------------------------------------------------------------------------------------------------
#remote server ime.usp.br backup sync
echo "Sending to $remote_host1:"
rsync -t -v -z -a --progress $output $user_remote_host1@$remote_host1:$remote_host1_dir
echo "Sent to $user_remote_host1@$remote_host1"
echo
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#remote server ybytu backup sync
echo "Sending to $remote_host2:"
ssh -t $user_remote_host1@$remote_host1 "rsync -t -v -z -a --progress $remote_host1_dir/$output $user_remote_host2@$remote_host2:$remote_host2_dir; rm -rf $output"
echo "Sent to $user_remote_host2@$remote_host2"
echo
#-------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------
# untar and compile
echo "Untar and compilation at $remote_host2"
ssh -t $user_remote_host1@$remote_host1 "ssh -t $user_remote_host2@$remote_host2 <<EOF
	cd $remote_host2_dir;
	tar -xvf $output;
	make clean;
	make;
	pwd;
	rm -rf $output;
EOF"
echo "Untar and compilation at $remote_host2 done."
#-------------------------------------------------------------------------------------------------------

# remove tar file
rm -rf $output
