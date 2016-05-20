PILOT=15
FINAL=65
U="2.3"
SIMS = 100

python $SCRIPTDIR/estimation_table.py $PILOT $FINAL $SIMS 'hcp' 'adaptive' $U $RESDIR $TABDIR
python $SCRIPTDIR/estimation_table.py $PILOT $FINAL $SIMS 'hcp' 'predictive' $U $RESDIR $TABDIR
python $SCRIPTDIR/estimation_table.py $PILOT $FINAL $SIMS 'sim' 'adaptive' $U $RESDIR $TABDIR
python $SCRIPTDIR/estimation_table.py $PILOT $FINAL $SIMS 'sim' 'predictive' $U $RESDIR $TABDIR

python estimation_powpre.py 15 80 250 47 'HCP' 'adaptive' 'u2' $RESDIR $TABDIR &
python estimation_powpre.py 15 80 250 47 'HCP' 'adaptive' 'u3' $RESDIR $TABDIR &
python estimation_powpre.py 15 80 250 47 'HCP' 'nonadaptive' 'u2' $RESDIR $TABDIR &
python estimation_powpre.py 15 80 250 47 'HCP' 'nonadaptive' 'u3' $RESDIR $TABDIR &
python estimation_powpre.py 15 61 250 16 'SIM' 'adaptive' 'u2' $RESDIR $TABDIR &
python estimation_powpre.py 15 61 250 16 'SIM' 'adaptive' 'u3' $RESDIR $TABDIR &
python estimation_powpre.py 15 61 250 16 'SIM' 'nonadaptive' 'u2' $RESDIR $TABDIR &
python estimation_powpre.py 15 61 250 16 'SIM' 'nonadaptive' 'u3' $RESDIR $TABDIR &

python estimation_powtrue.py 15 80 250 47 'HCP' 'adaptive' 'u2' $RESDIR $TABDIR &
python estimation_powtrue.py 15 80 250 47 'HCP' 'adaptive' 'u3' $RESDIR $TABDIR &
python estimation_powtrue.py 15 80 250 47 'HCP' 'nonadaptive' 'u2' $RESDIR $TABDIR &
python estimation_powtrue.py 15 80 250 47 'HCP' 'nonadaptive' 'u3' $RESDIR $TABDIR &
python estimation_powtrue.py 15 61 250 16 'SIM' 'adaptive' 'u2' $RESDIR $TABDIR &
python estimation_powtrue.py 15 61 250 16 'SIM' 'adaptive' 'u3' $RESDIR $TABDIR &
python estimation_powtrue.py 15 61 250 16 'SIM' 'nonadaptive' 'u2' $RESDIR $TABDIR &
python estimation_powtrue.py 15 61 250 16 'SIM' 'nonadaptive' 'u3' $RESDIR $TABDIR &

RESDIR="/scratch/users/jdurnez/"
TABDIR="/scratch/users/jdurnez/interim_results/"

SIMstart=3
SIMstop=4
python estimation_conditional_HCP.py 15 61 $SIMstart $SIMstop 16 'HCP' 'nonadaptive' 'u3' $RESDIR $TABDIR &

outfile="conditional_HCP_adaptive_u2.csv"
less conditional*_HCP_adaptive_u2.csv >> $outfile
header=$(head -n 1 $outfile)
sed s/$header//g $outfile > tmp && mv tmp $outfile
