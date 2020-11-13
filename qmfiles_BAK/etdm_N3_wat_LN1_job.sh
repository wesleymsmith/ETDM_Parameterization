
###START_PSI4     #DO NOT CHANGE
#$ -S /bin/bash  
#$ -cwd   
#$ -V     
#$ -pe smp 8   
#$ -l h_data=1200,h_rt=7200  
#$ -N ffp65  
#$ -o ffp65_out  
#$ -R y   
#$ -j y   

unset PSIDATADIR  
PSI4_ROOT=$HOME/bin/anaconda3/envs/psi4/bin/      
export PATH=$PSI4_ROOT/bin:$PATH  
export PYTHONPATH=$PSI4_ROOT/lib  
export PSI_SCRATCH=/tmp           

###########################
python etdm_N3_wat_LN1.py etdm_N3_wat_LN1.out               
###########################