cd ana/; ./p0.sh; cd ..
root -q -b -l 'macro_ana.C("14Oap")'
root -q -b -l 'macro_ana.C("14OCO2p")'
root -q -b -l 'macro_ana.C("14Oap_14N")'
root -q -b -l 'macro_ana.C("14OCO2p_14N")'

root macro_14O.C
