if [ $1 == 'gui' ]
then
    /home/cjrevelas/programs/comsol/multiphysics/bin/comsol
    exit
fi

if [ $1 == 'with_matlab' ]
then
    /home/cjrevelas/programs/comsol/multiphysics/bin/comsol server matlab
    exit
fi
