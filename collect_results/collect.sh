BASEFOLDER="/Users/Joke/Documents/Onderzoek/Studie_4_propow/InterimPower_Results"
OUTDIR=$(echo $BASEFOLDER\/tables/)

python estimation_table.py 15 80 250 'hcp' 'adaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_table.py 15 80 250 'hcp' 'adaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_table.py 15 80 250 'hcp' 'nonadaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_table.py 15 80 250 'hcp' 'nonadaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_table.py 15 61 250 'sim' 'adaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_table.py 15 61 250 'sim' 'adaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_table.py 15 61 250 'sim' 'nonadaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_table.py 15 61 250 'sim' 'nonadaptive' 'u3' $BASEFOLDER $OUTDIR &

python estimation_powpre.py 15 80 250 47 'hcp' 'adaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_powpre.py 15 80 250 47 'hcp' 'adaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_powpre.py 15 80 250 47 'hcp' 'nonadaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_powpre.py 15 80 250 47 'hcp' 'nonadaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_powpre.py 15 61 250 16 'sim' 'adaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_powpre.py 15 61 250 16 'sim' 'adaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_powpre.py 15 61 250 16 'sim' 'nonadaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_powpre.py 15 61 250 16 'sim' 'nonadaptive' 'u3' $BASEFOLDER $OUTDIR &

python estimation_powtrue.py 15 80 250 47 'hcp' 'adaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_powtrue.py 15 80 250 47 'hcp' 'adaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_powtrue.py 15 80 250 47 'hcp' 'nonadaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_powtrue.py 15 80 250 47 'hcp' 'nonadaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_powtrue.py 15 61 250 16 'sim' 'adaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_powtrue.py 15 61 250 16 'sim' 'adaptive' 'u3' $BASEFOLDER $OUTDIR &
python estimation_powtrue.py 15 61 250 16 'sim' 'nonadaptive' 'u2' $BASEFOLDER $OUTDIR &
python estimation_powtrue.py 15 61 250 16 'sim' 'nonadaptive' 'u3' $BASEFOLDER $OUTDIR &

BASEFOLDER="/scratch/users/jdurnez/"
OUTDIR="/scratch/users/jdurnez/interim_results/"

simstart=3
simstop=4
python estimation_conditional_hcp.py 15 61 $simstart $simstop 16 'hcp' 'nonadaptive' 'u3' $BASEFOLDER $OUTDIR &

outfile="conditional_hcp_adaptive_u2.csv"
less conditional*_hcp_adaptive_u2.csv >> $outfile
header=$(head -n 1 $outfile)
sed s/$header//g $outfile > tmp && mv tmp $outfile
