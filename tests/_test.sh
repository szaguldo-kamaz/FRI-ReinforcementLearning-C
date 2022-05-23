mkdir -p tmp
cd tmp

#command -v mpirun >/dev/null 2>&1 && RUN="mpirun -n 4"

#time $RUN $APP

time $APP -q

diff $RB_NEW $RB_ORIG

if [ $? -eq 0 ]; then
        echo Valid
else
        echo Invalid
        vimdiff $RB_NEW $RB_ORIG
fi
