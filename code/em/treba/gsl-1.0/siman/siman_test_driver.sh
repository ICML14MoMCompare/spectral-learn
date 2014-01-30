#! /bin/sh

# assume good result from tests; increment it if any test fails
EXIT_STATUS=0

echo -n "running simulated annealing test #1 (travelling salesman problem)"
./siman_test > siman_test.out1
SECOND_LAST_ENERGY=`tail -2 siman_test.out1 | head -1 | awk '{print $4}'`
LAST_ENERGY=`tail -1 siman_test.out1 | awk '{print $4}'`
# echo " " $SECOND_LAST_ENERGY $LAST_ENERGY
if [ $SECOND_LAST_ENERGY = $LAST_ENERGY ];
then
    echo "... converged -- PASS"
else
    echo "... did NOT converge -- FAIL"
    EXIT_STATUS=`expr $EXIT_STATUS + 1`
    echo $EXIT_STATUS
fi

echo -n "running simulated annealing test #2 (travelling salesman problem)"
GSL_RNG_SEED=12345 ./siman_test > siman_test.out2
SECOND_LAST_ENERGY=`tail -2 siman_test.out2 | head -1 | awk '{print $4}'`
LAST_ENERGY=`tail -1 siman_test.out2 | awk '{print $4}'`
# echo " " $SECOND_LAST_ENERGY $LAST_ENERGY
if [ $SECOND_LAST_ENERGY = $LAST_ENERGY ];
then
    echo "... converged -- PASS"
else
    echo "... did NOT converge -- FAIL"
    EXIT_STATUS=`expr $EXIT_STATUS + 1`
    echo $EXIT_STATUS
fi

exit $EXIT_STATUS
