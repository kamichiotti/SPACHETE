#RB I don't think this is appending to FJ, just overwriting it

for file in *
do
FJ_input=${file}
done

printf "$FJ_input\n"
