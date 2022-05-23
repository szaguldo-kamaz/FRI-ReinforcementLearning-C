You can find the original and valid outputs of the demos in the "orig" directory.
Execute the test_*.sh bash scripts to check it against your results.

Example:

cd <build_dir>/tests
./test_acrobot.sh

To manually check the correctness, run the demo, and diff the output:

cd <build_dir>/examples/acrobot
time ./acrobot && diff ./acrobot.frirlrb.txt <this_dir>/original/frirl_example_cartpole.frirlrb.txt && echo $?
