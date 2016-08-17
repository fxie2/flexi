#!/bin/bash

testcases=(  '1.1'    '1.2'      '1.3'        '2.1'     '2.2'      '2.3'      '3.1'    '3.3'   '4.1'    '5.1'    '5.2'     '5.3'      '6.1'    '7.1'       '7.2'     '7.3')
length=(     '1.e-5'  '1.e-4'    '1.e-3'      '1.e-6'   '1.e-4'    '1.e-3'    '6.8e-7' '1.e-3' '3.4e-7' '2.e-7'  '1.e-4'   '1.e-3'    '1.7e-7' '1.3075e-7' '1.e-4'   '1.e-3')
pressure=(   '101325' '20000'    '2000'       '101325'  '2000'     '400'      '101325' '400'   '101325' '101325' '400'     '100'      '101325' '101325'    '400'     '100')
temperature=('300'    '600'      '600'        '300'     '600'      '1200'     '300'    '1740'  '300'    '300'    '600'     '1500'     '300'    '300'       '900'     '2250')
density=(    '1.225'  '0.1165'   '0.01165'    '1.225'   '0.01165'  '0.001165' '1.225'  '0.0008' '1.225' '1.225'  '0.00233' '0.000235' '1.225'  '1.225'     '0.00155' '0.000158')
analyze=(    '5.e-6'  '5.e-5'    '5.e-4'      '2.e-7'   '2.e-5'    '2.e-4'    '1.e-7'  '1.e-4' '2.e-8'  '1.e-8'  '5.e-6'   '5.e-5'    '1.e-8'  '5.e-9'     '2.e-6'   '2.e-5')

ntests=${#testcases[@]}

for ((icase=0;icase<${#testcases[@]};++icase)); do
  mycase=${testcases[$icase]}
  mylen=${length[$icase]}
  mypresout=${pressure[$icase]}
  mypresin=`echo "scale=10; $mypresout + 5.0" | bc`
  mytemp=${temperature[$icase]}
  mydensout=${density[$icase]}
  mydensin=`echo "scale=10; $mypresin / $mypresout * $mydensout" | bc`  # isothermal
  myanadt=${analyze[$icase]}
  myR=`echo "scale=10; $mypresout/($mydensout*$mytemp)" | bc`
  myvisc=`echo "scale=10; 5.0/16.0*sqrt(3.14159265359/($myR*$mytemp))*0.001380648813*$mytemp/(3.14159265359*3.7*3.7)" | bc` 

  echo "-------------------------------------------------------------------------"
  echo "Computing case number $mycase"
  echo "-------------------------------------------------------------------------"
  rm -rf CASE_${mycase}
  mkdir CASE_${mycase}
  cd CASE_${mycase}
  cat ../preproc_sample.ini   |  sed -e 's/FIGICASE/'$mycase'/g'       | \
                                 sed -e 's/FIGILEN/'$mylen'/g'           \
                                 >> preproc_${mycase}.ini
  ../../../bin/preproctool preproc_${mycase}.ini | tee CASE_${mycase}.out

  cat ../parameter_sample.ini |  sed -e 's/FIGICASE/'$mycase'/g'       | \
                                 sed -e 's/FIGILEN/'$mylen'/g'         | \
                                 sed -e 's/FIGIPRESIN/'$mypresin'/g'   | \
                                 sed -e 's/FIGIPRESOUT/'$mypresout'/g' | \
                                 sed -e 's/FIGIDENSOUT/'$mydensout'/g' | \
                                 sed -e 's/FIGIDENSIN/'$mydensin'/g'   | \
                                 sed -e 's/FIGIR/'$myR'/g'             | \
                                 sed -e 's/FIGIVISC/'$myvisc'/g'       | \
                                 sed -e 's/FIGIANADT/'$myanadt'/g'       \
                                 >> parameter_${mycase}.ini

  #timeout --signal=KILL 2000s mpirun -np 5 ../../../bin/flexi parameter_${mycase}.ini | tee CASE_${mycase}.out
  timeout --signal=KILL 3600s mpirun -np 5 ../../../bin/flexi parameter_${mycase}.ini | tee CASE_${mycase}.out

  cd ../

done
