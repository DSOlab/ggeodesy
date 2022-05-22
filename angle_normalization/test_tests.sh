#! /bin/bash

progs=(nang10.out nang20.out nang30.out nang40.out nang12.out nang22.out nang32.out nang42.out)
pave=(0 0 0 0 0 0 0 0)
pn=(0 0 0 0 0 0 0 0)

echo "performing tests ..."
for i in `seq 1 1000` ; do
  tnr=$(python -c "import random; print(random.randint(0,7))")
  #echo "performing test for algorithm $tnr, average=${pave[$tnr]}, N=${pn[$tnr]}"
  tm=$(./${progs[$tnr]} | sed -rn 's/.*Running time: ([0-9]*) .*/\1/p')
  nv=$(python -c "print(${pave[$tnr]} + (-${pave[$tnr]}+$tm)/(${pn[$tnr]}+1e0))")
  pave[$tnr]=$nv
  pn[$tnr]=$(python -c "print(${pn[$tnr]}+1)")
  #echo "new values for algorithm $tnr, ${pave[$tnr]} and ${pn[$tnr]}"
  #echo "${pave[*]} // ${pn[*]}"
done

for i in `seq 0 8` ; do
  echo "Program ${progs[$i]} tests: ${pn[$i]} Average running time: ${pave[$i]}"
done
